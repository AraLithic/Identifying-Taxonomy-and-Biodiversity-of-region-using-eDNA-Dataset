"""
TerraTrace backend by Umair
Deterministic analysis of processed eDNA / metabarcoding CSV datasets.

No AI. No external API calls. No database. Everything is computed in
memory, directly from the uploaded CSV, using pandas / numpy / scipy.
"""

import io
import math
import random
from datetime import datetime, timedelta
from typing import Optional

import numpy as np
import pandas as pd
from fastapi import FastAPI, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, HTMLResponse
from scipy.spatial.distance import pdist, squareform

# App setup

app = FastAPI(title="TerraTrace", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

MAX_UPLOAD_BYTES = 25 * 1024 * 1024  # 25 MB se bada nahi kuch
MAX_PAIRWISE_SITES = 60  # protect pairwise matrix computation from blowing up shit out of my PC

# Column detection

TAXON_CANDIDATES = ["species", "scientific_name", "taxon", "organism"]
SAMPLE_CANDIDATES = ["sample_id", "sample", "sample_name"]
SITE_CANDIDATES = ["site", "location", "sampling_site", "station"]
ABUNDANCE_CANDIDATES = ["reads", "read_count", "abundance", "count"]
SEQUENCE_CANDIDATES = ["sequence", "seq", "dna_sequence"]
LAT_CANDIDATES = ["latitude", "lat"]
LON_CANDIDATES = ["longitude", "lon", "lng", "long"]
DATE_CANDIDATES = ["date", "sampling_date", "collection_date"]

RANK_COLUMNS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def _find_column(columns_lower_map, candidates):
    """Return the original column name matching the first candidate found."""
    for candidate in candidates:
        if candidate in columns_lower_map:
            return columns_lower_map[candidate]
    return None


def detect_columns(df: pd.DataFrame) -> dict:
    """Conservatively detect the semantic role of each column."""
    columns_lower_map = {c.lower().strip(): c for c in df.columns}

    detected = {
        "taxon": _find_column(columns_lower_map, TAXON_CANDIDATES),
        "sample": _find_column(columns_lower_map, SAMPLE_CANDIDATES),
        "site": _find_column(columns_lower_map, SITE_CANDIDATES),
        "abundance": _find_column(columns_lower_map, ABUNDANCE_CANDIDATES),
        "sequence": _find_column(columns_lower_map, SEQUENCE_CANDIDATES),
        "latitude": _find_column(columns_lower_map, LAT_CANDIDATES),
        "longitude": _find_column(columns_lower_map, LON_CANDIDATES),
        "date": _find_column(columns_lower_map, DATE_CANDIDATES),
    }

    rank_columns = {}
    for rank in RANK_COLUMNS:
        if rank in columns_lower_map:
            rank_columns[rank] = columns_lower_map[rank]
    detected["ranks"] = rank_columns

    # If no dedicated "taxon" column was found, but a "species" rank column
    # exists, it can double as the taxon identifier.
    if detected["taxon"] is None and "species" in rank_columns:
        detected["taxon"] = rank_columns["species"]

    detected["confident"] = detected["taxon"] is not None
    return detected

# Data quality

def assess_quality(df: pd.DataFrame, columns: dict) -> dict:
    warnings = []
    n_rows, n_cols = df.shape

    missing_values = int(df.isna().sum().sum())
    duplicate_rows = int(df.duplicated().sum())

    invalid_abundance = 0
    if columns["abundance"]:
        col = df[columns["abundance"]]
        numeric = pd.to_numeric(col, errors="coerce")
        invalid_abundance = int(numeric.isna().sum() - col.isna().sum())
        invalid_abundance += int((numeric < 0).sum())
        if invalid_abundance > 0:
            warnings.append(
                f"{invalid_abundance} row(s) had a non-numeric or negative "
                f"value in the abundance column '{columns['abundance']}'."
            )

    invalid_coordinates = 0
    if columns["latitude"] and columns["longitude"]:
        lat = pd.to_numeric(df[columns["latitude"]], errors="coerce")
        lon = pd.to_numeric(df[columns["longitude"]], errors="coerce")
        bad_lat = (lat.isna()) | (lat < -90) | (lat > 90)
        bad_lon = (lon.isna()) | (lon < -180) | (lon > 180)
        invalid_coordinates = int((bad_lat | bad_lon).sum())
        if invalid_coordinates > 0:
            warnings.append(
                f"{invalid_coordinates} row(s) had missing or out-of-range "
                f"latitude/longitude values."
            )

    invalid_sequences = 0
    if columns["sequence"]:
        valid_chars = set("ACGTN")
        seqs = df[columns["sequence"]].astype(str)

        def is_invalid(s):
            s_upper = s.strip().upper()
            if s_upper in ("", "NAN", "NONE"):
                return True
            return any(ch not in valid_chars for ch in s_upper)

        invalid_sequences = int(seqs.apply(is_invalid).sum())
        if invalid_sequences > 0:
            warnings.append(
                f"{invalid_sequences} row(s) had a sequence with characters "
                f"outside A/C/G/T/N or an empty sequence."
            )

    if duplicate_rows > 0:
        warnings.append(
            f"{duplicate_rows} fully duplicate row(s) were detected and kept "
            f"in the dataset (not automatically removed)."
        )

    if missing_values > 0:
        warnings.append(f"{missing_values} missing value(s) were found across the dataset.")

    if not columns["taxon"]:
        warnings.append(
            "No taxon/species column could be confidently identified. "
            "Please map a column manually."
        )

    return {
        "rows": int(n_rows),
        "columns": int(n_cols),
        "missing_values": missing_values,
        "duplicate_rows": duplicate_rows,
        "invalid_abundance_values": invalid_abundance,
        "invalid_coordinates": invalid_coordinates,
        "invalid_sequences": invalid_sequences,
        "warnings": warnings,
    }


# Taxonomy (NCERT padho agar samjh na aye)

def analyze_taxonomy(df: pd.DataFrame, columns: dict) -> dict:
    if not columns["taxon"]:
        return {"available": False, "reason": "No taxon/species column identified."}

    taxon_col = columns["taxon"]
    taxa_series = df[taxon_col].dropna().astype(str).str.strip()
    taxa_series = taxa_series[taxa_series != ""]

    unique_taxa = sorted(taxa_series.unique().tolist())
    richness = len(unique_taxa)

    # Detection/abundance counts per taxon
    if columns["abundance"]:
        work = df[[taxon_col, columns["abundance"]]].copy()
        work[columns["abundance"]] = pd.to_numeric(work[columns["abundance"]], errors="coerce")
        work = work.dropna()
        work = work[work[columns["abundance"]] >= 0]
        counts = work.groupby(taxon_col)[columns["abundance"]].sum().sort_values(ascending=False)
        count_basis = "abundance"
    else:
        counts = taxa_series.value_counts()
        count_basis = "presence_count"

    taxon_table = []
    ranks_present = list(columns["ranks"].keys())
    for taxon_name, count in counts.items():
        row_match = df[df[taxon_col].astype(str).str.strip() == taxon_name]
        entry = {"taxon": taxon_name, "count": float(count), "count_basis": count_basis}
        for rank in RANK_COLUMNS:
            rank_col = columns["ranks"].get(rank)
            if rank_col and not row_match.empty:
                val = row_match.iloc[0][rank_col]
                entry[rank] = "Not provided" if pd.isna(val) or str(val).strip() == "" else str(val).strip()
            else:
                entry[rank] = "Not provided"
        taxon_table.append(entry)

    return {
        "available": True,
        "richness": richness,
        "ranks_present": ranks_present,
        "count_basis": count_basis,
        "taxon_table": taxon_table,
    }

# Biodiversity metrics mat samjho

def shannon_index(proportions: np.ndarray) -> float:
    p = proportions[proportions > 0]
    return float(-np.sum(p * np.log(p)))


def simpson_diversity(proportions: np.ndarray) -> float:
    p = proportions[proportions > 0]
    d = float(np.sum(p ** 2))
    return 1.0 - d


def pielou_evenness(shannon: float, richness: int) -> Optional[float]:
    if richness > 1:
        return float(shannon / math.log(richness))
    return None


def margalef_richness(richness: int, total_individuals: float) -> Optional[float]:
    if richness > 0 and total_individuals > 1:
        return float((richness - 1) / math.log(total_individuals))
    return None


def chao1_estimate(counts: np.ndarray) -> Optional[float]:
    """Chao1 richness estimator. Requires integer-like abundance counts."""
    s_obs = int(np.sum(counts > 0))
    f1 = int(np.sum(counts == 1))  # singletons
    f2 = int(np.sum(counts == 2))  # doubletons
    if s_obs == 0:
        return None
    if f2 > 0:
        return float(s_obs + (f1 ** 2) / (2 * f2))
    if f1 > 0:
        # bias-corrected form when there are no doubletons
        return float(s_obs + (f1 * (f1 - 1)) / 2)
    return float(s_obs)


def compute_biodiversity(taxon_counts: pd.Series, abundance_mode: bool) -> dict:
    counts = taxon_counts.values.astype(float)
    richness = int(np.sum(counts > 0))
    total = float(np.sum(counts))

    if total <= 0 or richness == 0:
        return {"available": False, "reason": "No taxon counts available for this group."}

    proportions = counts / total
    shannon = shannon_index(proportions)
    simpson = simpson_diversity(proportions)
    evenness = pielou_evenness(shannon, richness)
    margalef = margalef_richness(richness, total)

    result = {
        "available": True,
        "observed_richness": richness,
        "shannon_diversity": round(shannon, 4),
        "simpson_diversity": round(simpson, 4),
        "pielou_evenness": round(evenness, 4) if evenness is not None else None,
        "margalef_richness": round(margalef, 4) if margalef is not None else None,
        "total_count": total,
        "count_basis": "abundance" if abundance_mode else "presence_count",
    }

    if abundance_mode:
        # Chao1 only makes sense on integer-like abundance/read counts
        is_integer_like = np.allclose(counts, np.round(counts))
        if is_integer_like:
            chao1 = chao1_estimate(np.round(counts))
            result["chao1_estimated_richness"] = round(chao1, 4) if chao1 is not None else None
        else:
            result["chao1_estimated_richness"] = None
            result["chao1_reason"] = "Abundance values are not integer read/detection counts."
    else:
        result["chao1_estimated_richness"] = None
        result["chao1_reason"] = "Chao1 requires abundance/read-count data, which was not available."

    return result

# Community similarity

def build_site_taxon_matrix(df: pd.DataFrame, columns: dict, abundance_mode: bool):
    taxon_col = columns["taxon"]
    site_col = columns["site"]

    work = df[[site_col, taxon_col]].copy()
    if abundance_mode:
        work["value"] = pd.to_numeric(df[columns["abundance"]], errors="coerce")
    else:
        work["value"] = 1.0

    work = work.dropna(subset=[site_col, taxon_col])
    work["value"] = work["value"].fillna(0)
    work = work[work["value"] >= 0]

    matrix = work.pivot_table(
        index=site_col, columns=taxon_col, values="value", aggfunc="sum", fill_value=0
    )
    return matrix


def compute_community_similarity(df: pd.DataFrame, columns: dict, abundance_mode: bool) -> dict:
    if not columns["site"] or not columns["taxon"]:
        return {"available": False, "reason": "A site column and a taxon column are both required."}

    matrix = build_site_taxon_matrix(df, columns, abundance_mode)
    n_sites = matrix.shape[0]

    if n_sites < 2:
        return {"available": False, "reason": "At least two sites are required to compare community composition."}

    if n_sites > MAX_PAIRWISE_SITES:
        return {
            "available": False,
            "reason": f"Too many sites ({n_sites}) for a pairwise comparison in this prototype "
                      f"(limit {MAX_PAIRWISE_SITES}). Please analyze a smaller dataset.",
        }

    sites = matrix.index.tolist()
    presence_matrix = (matrix > 0).astype(int)

    result = {"available": True, "sites": sites}

    # Jaccard and Sorensen use presence/absence, always available if a
    # taxon+site column exist.
    jaccard_condensed = pdist(presence_matrix.values, metric="jaccard")
    jaccard_matrix = squareform(jaccard_condensed)
    result["jaccard"] = jaccard_matrix.round(4).tolist()

    # Sorensen-Dice: 2 * shared / (a + b)
    sorensen_condensed = pdist(presence_matrix.values, metric="dice")
    sorensen_matrix = squareform(sorensen_condensed)
    result["sorensen"] = sorensen_matrix.round(4).tolist()

    if abundance_mode:
        bray_condensed = pdist(matrix.values, metric="braycurtis")
        bray_matrix = np.nan_to_num(squareform(bray_condensed), nan=0.0)
        result["bray_curtis"] = bray_matrix.round(4).tolist()
    else:
        result["bray_curtis"] = None
        result["bray_curtis_reason"] = "Bray-Curtis requires abundance/read-count data, which was not available."

    return result

# Site-level analysis

def analyze_sites(df: pd.DataFrame, columns: dict, abundance_mode: bool) -> dict:
    if not columns["site"]:
        return {"available": False, "reason": "No site/location column identified."}
    if not columns["taxon"]:
        return {"available": False, "reason": "No taxon column identified."}

    site_col = columns["site"]
    taxon_col = columns["taxon"]

    sites_out = []
    for site_name, group in df.groupby(site_col):
        n_samples = int(group[columns["sample"]].nunique()) if columns["sample"] else int(len(group))
        taxa = group[taxon_col].dropna().astype(str)
        taxa = taxa[taxa.str.strip() != ""]

        if abundance_mode:
            g = group.copy()
            g[columns["abundance"]] = pd.to_numeric(g[columns["abundance"]], errors="coerce")
            g = g.dropna(subset=[columns["abundance"]])
            counts = g.groupby(taxon_col)[columns["abundance"]].sum()
            total = float(counts.sum())
        else:
            counts = taxa.value_counts()
            total = float(counts.sum())

        bio = compute_biodiversity(counts, abundance_mode) if len(counts) > 0 else {"available": False}

        sites_out.append({
            "site": str(site_name),
            "n_samples": n_samples,
            "richness": int((counts > 0).sum()) if len(counts) else 0,
            "shannon_diversity": bio.get("shannon_diversity") if bio.get("available") else None,
            "simpson_diversity": bio.get("simpson_diversity") if bio.get("available") else None,
            "total_count": total,
            "count_basis": "abundance" if abundance_mode else "presence_count",
        })

    sites_out.sort(key=lambda x: x["site"])
    return {"available": True, "sites": sites_out}

# Sequence QC (optional)

def analyze_sequence_qc(df: pd.DataFrame, columns: dict) -> dict:
    if not columns["sequence"]:
        return {"available": False, "reason": "No DNA sequence column was detected."}

    seq_col = columns["sequence"]
    raw_seqs = df[seq_col].dropna().astype(str).str.strip()
    raw_seqs = raw_seqs[raw_seqs != ""]

    if len(raw_seqs) == 0:
        return {"available": False, "reason": "The sequence column contained no usable values."}

    valid_chars = set("ACGTN")
    lengths = []
    gc_percentages = []
    invalid_count = 0
    n_count_total = 0

    for s in raw_seqs:
        s_upper = s.upper()
        if any(ch not in valid_chars for ch in s_upper):
            invalid_count += 1
            continue
        length = len(s_upper)
        lengths.append(length)
        g = s_upper.count("G")
        c = s_upper.count("C")
        n_count_total += s_upper.count("N")
        gc_percentages.append((g + c) / length * 100 if length > 0 else 0)

    duplicate_count = int(raw_seqs.duplicated().sum())

    if not lengths:
        return {
            "available": True,
            "valid_sequences": 0,
            "invalid_sequences": invalid_count,
            "duplicate_sequences": duplicate_count,
            "length_stats": None,
            "gc_stats": None,
            "ambiguous_base_count": n_count_total,
            "length_distribution": [],
            "gc_distribution": [],
        }

    lengths_arr = np.array(lengths)
    gc_arr = np.array(gc_percentages)

    return {
        "available": True,
        "valid_sequences": len(lengths),
        "invalid_sequences": invalid_count,
        "duplicate_sequences": duplicate_count,
        "ambiguous_base_count": n_count_total,
        "length_stats": {
            "min": int(lengths_arr.min()),
            "max": int(lengths_arr.max()),
            "mean": round(float(lengths_arr.mean()), 2),
            "median": round(float(np.median(lengths_arr)), 2),
        },
        "gc_stats": {
            "min": round(float(gc_arr.min()), 2),
            "max": round(float(gc_arr.max()), 2),
            "mean": round(float(gc_arr.mean()), 2),
            "median": round(float(np.median(gc_arr)), 2),
        },
        "length_distribution": lengths_arr.tolist(),
        "gc_distribution": np.round(gc_arr, 2).tolist(),
    }

# Spatial analysis (optional) Gand faad Math

def analyze_spatial(df: pd.DataFrame, columns: dict) -> dict:
    if not columns["latitude"] or not columns["longitude"]:
        return {"available": False, "reason": "No valid latitude/longitude columns were detected."}

    lat = pd.to_numeric(df[columns["latitude"]], errors="coerce")
    lon = pd.to_numeric(df[columns["longitude"]], errors="coerce")
    valid = (~lat.isna()) & (~lon.isna()) & lat.between(-90, 90) & lon.between(-180, 180)

    if valid.sum() == 0:
        return {"available": False, "reason": "No valid latitude/longitude coordinate pairs were found."}

    points = []
    site_col = columns["site"]
    taxon_col = columns["taxon"]
    subset = df[valid]
    for _, row in subset.iterrows():
        point = {"latitude": float(row[columns["latitude"]]), "longitude": float(row[columns["longitude"]])}
        if site_col:
            point["site"] = str(row[site_col]) if not pd.isna(row[site_col]) else None
        if taxon_col:
            point["taxon"] = str(row[taxon_col]) if not pd.isna(row[taxon_col]) else None
        points.append(point)

    return {"available": True, "points": points, "n_valid": int(valid.sum()), "n_invalid": int((~valid).sum())}

# Temporal analysis (optional) Temporal Lobe jaisa hai function

def analyze_temporal(df: pd.DataFrame, columns: dict) -> dict:
    if not columns["date"]:
        return {"available": False, "reason": "No valid sampling-date column was detected."}

    dates = pd.to_datetime(df[columns["date"]], errors="coerce")
    valid = ~dates.isna()

    if valid.sum() == 0:
        return {"available": False, "reason": "No valid dates could be parsed from the date column."}

    taxon_col = columns["taxon"]
    work = df[valid].copy()
    work["_parsed_date"] = dates[valid].dt.date.astype(str)

    samples_over_time = work.groupby("_parsed_date").size().sort_index()
    series = [{"date": d, "sample_count": int(c)} for d, c in samples_over_time.items()]

    taxa_over_time = None
    if taxon_col:
        taxa_grouped = work.groupby("_parsed_date")[taxon_col].nunique().sort_index()
        taxa_over_time = [{"date": d, "unique_taxa": int(c)} for d, c in taxa_grouped.items()]

    return {
        "available": True,
        "samples_over_time": series,
        "taxa_over_time": taxa_over_time,
        "n_valid_dates": int(valid.sum()),
        "n_invalid_dates": int((~valid).sum()),
    }

# Core analysis pipeline(Danger: Don't touch)

def run_analysis(df: pd.DataFrame, is_demo: bool = False) -> dict:
    warnings = []

    if df.empty:
        return {"error": "The uploaded CSV contains no rows of data."}

    columns = detect_columns(df)

    if not columns["taxon"]:
        warnings.append(
            "TerraTrace could not confidently identify a taxon column. "
            "Please use the column mapping panel to select one manually."
        )

    quality = assess_quality(df, columns)

    abundance_mode = columns["abundance"] is not None
    if not abundance_mode:
        warnings.append("Abundance information was not available. Presence/absence analysis is being used.")

    taxonomy = analyze_taxonomy(df, columns)

    biodiversity = {"available": False, "reason": "No taxon column identified."}
    if columns["taxon"]:
        taxon_col = columns["taxon"]
        if abundance_mode:
            work = df[[taxon_col, columns["abundance"]]].copy()
            work[columns["abundance"]] = pd.to_numeric(work[columns["abundance"]], errors="coerce")
            work = work.dropna()
            work = work[work[columns["abundance"]] >= 0]
            counts = work.groupby(taxon_col)[columns["abundance"]].sum()
        else:
            taxa_series = df[taxon_col].dropna().astype(str).str.strip()
            taxa_series = taxa_series[taxa_series != ""]
            counts = taxa_series.value_counts()
        biodiversity = compute_biodiversity(counts, abundance_mode)

    sites = analyze_sites(df, columns, abundance_mode)
    community = compute_community_similarity(df, columns, abundance_mode)
    sequence_qc = analyze_sequence_qc(df, columns)
    spatial = analyze_spatial(df, columns)
    temporal = analyze_temporal(df, columns)

    detected_columns_summary = {
        "taxon": columns["taxon"],
        "sample": columns["sample"],
        "site": columns["site"],
        "abundance": columns["abundance"],
        "sequence": columns["sequence"],
        "latitude": columns["latitude"],
        "longitude": columns["longitude"],
        "date": columns["date"],
        "ranks": columns["ranks"],
        "all_columns": df.columns.tolist(),
        "confident": columns["confident"],
        "analysis_mode": "abundance" if abundance_mode else "presence_absence",
    }

    all_warnings = quality["warnings"] + warnings

    return {
        "is_demo": is_demo,
        "dataset": {
            "rows": int(df.shape[0]),
            "columns": int(df.shape[1]),
        },
        "columns": detected_columns_summary,
        "quality": quality,
        "taxonomy": taxonomy,
        "biodiversity": biodiversity,
        "sites": sites,
        "community": community,
        "sequence_qc": sequence_qc,
        "spatial": spatial,
        "temporal": temporal,
        "warnings": all_warnings,
    }

# Demo dataset

def generate_demo_dataframe() -> pd.DataFrame:
    rng = random.Random(42)

    taxa_pool = [
        ("Animalia", "Chordata", "Actinopterygii", "Cypriniformes", "Cyprinidae", "Labeo", "Labeo rohita"),
        ("Animalia", "Chordata", "Actinopterygii", "Cypriniformes", "Cyprinidae", "Catla", "Catla catla"),
        ("Animalia", "Chordata", "Actinopterygii", "Siluriformes", "Bagridae", "Mystus", "Mystus tengara"),
        ("Animalia", "Mollusca", "Bivalvia", "Unionida", "Unionidae", "Lamellidens", "Lamellidens marginalis"),
        ("Animalia", "Arthropoda", "Malacostraca", "Decapoda", "Palaemonidae", "Macrobrachium", "Macrobrachium lamarrei"),
        ("Animalia", "Chordata", "Amphibia", "Anura", "Dicroglossidae", "Euphlyctis", "Euphlyctis cyanophlyctis"),
        ("Plantae", "Charophyta", "Chlorophyceae", "Charales", "Characeae", "Chara", "Chara vulgaris"),
        ("Animalia", "Rotifera", "Eurotatoria", "Ploima", "Brachionidae", "Brachionus", "Brachionus calyciflorus"),
        ("Chromista", "Bacillariophyta", "Bacillariophyceae", "Naviculales", "Naviculaceae", "Navicula", "Navicula radiosa"),
        ("Animalia", "Chordata", "Actinopterygii", "Perciformes", "Nandidae", "Nandus", "Nandus nandus"),
        ("Animalia", "Arthropoda", "Insecta", "Odonata", "Libellulidae", "Orthetrum", "Orthetrum sabina"),
        ("Animalia", "Annelida", "Clitellata", "Haplotaxida", "Naididae", "Tubifex", "Tubifex tubifex"),
    ]

    sites = ["Site-A", "Site-B", "Site-C", "Site-D"]
    site_coords = {
        "Site-A": (20.2961, 85.8245),
        "Site-B": (20.3560, 85.8830),
        "Site-C": (20.2500, 85.7900),
        "Site-D": (20.4100, 85.9200),
    }

    bases = "ACGT"
    start_date = datetime(2025, 1, 6)

    rows = []
    sample_counter = 1
    for site in sites:
        base_lat, base_lon = site_coords[site]
        n_samples = rng.randint(3, 5)
        for s in range(n_samples):
            sample_id = f"S{sample_counter:03d}"
            sample_counter += 1
            sample_date = start_date + timedelta(days=rng.randint(0, 90))
            n_taxa_in_sample = rng.randint(4, 8)
            chosen = rng.sample(taxa_pool, n_taxa_in_sample)
            for taxon in chosen:
                kingdom, phylum, cls, order, family, genus, species = taxon
                reads = rng.randint(5, 500)
                seq_len = rng.randint(140, 220)
                sequence = "".join(rng.choice(bases) for _ in range(seq_len))
                rows.append({
                    "sample_id": sample_id,
                    "site": site,
                    "species": species,
                    "kingdom": kingdom,
                    "phylum": phylum,
                    "class": cls,
                    "order": order,
                    "family": family,
                    "genus": genus,
                    "reads": reads,
                    "sequence": sequence,
                    "latitude": round(base_lat + rng.uniform(-0.01, 0.01), 6),
                    "longitude": round(base_lon + rng.uniform(-0.01, 0.01), 6),
                    "date": sample_date.strftime("%Y-%m-%d"),
                })

    return pd.DataFrame(rows)

# API endpoints(Chuna mat sab tu jayega, bhut marunga)

@app.get("/")
def root():
    return HTMLResponse(
        "<html><body><h3>TerraTrace backend is running.</h3>"
        "<p>Open frontend/index.html in your browser to use the application.</p>"
        "</body></html>"
    )


@app.get("/api/health")
def health():
    return {"status": "ok", "service": "TerraTrace"}


@app.get("/api/demo")
def demo():
    try:
        df = generate_demo_dataframe()
        result = run_analysis(df, is_demo=True)
        return JSONResponse(content=result)
    except Exception:
        return JSONResponse(
            status_code=500,
            content={"error": "The demo dataset could not be analyzed. Please try again."},
        )


@app.post("/api/analyze")
async def analyze(file: UploadFile = File(...)):
    filename = file.filename or ""
    if not filename.lower().endswith(".csv"):
        return JSONResponse(status_code=400, content={"error": "Please upload a CSV file."})

    raw = await file.read()

    if len(raw) == 0:
        return JSONResponse(status_code=400, content={"error": "The uploaded file is empty."})

    if len(raw) > MAX_UPLOAD_BYTES:
        return JSONResponse(
            status_code=400,
            content={"error": f"File is too large. Maximum upload size is {MAX_UPLOAD_BYTES // (1024*1024)} MB."},
        )

    try:
        text_buffer = io.StringIO(raw.decode("utf-8", errors="replace"))
        df = pd.read_csv(text_buffer)
    except Exception:
        return JSONResponse(
            status_code=400,
            content={"error": "The file could not be parsed as a valid CSV. Please check its formatting."},
        )

    if df.empty or df.shape[1] == 0:
        return JSONResponse(status_code=400, content={"error": "The uploaded CSV contains no usable data."})

    try:
        result = run_analysis(df, is_demo=False)
    except Exception:
        return JSONResponse(
            status_code=500,
            content={"error": "An unexpected error occurred while analyzing the dataset. Please check the file and try again."},
        )

    return JSONResponse(content=result)


@app.post("/api/analyze-mapped")
async def analyze_mapped(
    file: UploadFile = File(...),
    taxon_column: Optional[str] = None,
    sample_column: Optional[str] = None,
    site_column: Optional[str] = None,
    abundance_column: Optional[str] = None,
):
    """Re-run analysis using a manual column mapping supplied by the user."""
    filename = file.filename or ""
    if not filename.lower().endswith(".csv"):
        return JSONResponse(status_code=400, content={"error": "Please upload a CSV file."})

    raw = await file.read()
    if len(raw) == 0:
        return JSONResponse(status_code=400, content={"error": "The uploaded file is empty."})

    try:
        text_buffer = io.StringIO(raw.decode("utf-8", errors="replace"))
        df = pd.read_csv(text_buffer)
    except Exception:
        return JSONResponse(status_code=400, content={"error": "The file could not be parsed as a valid CSV."})

    if df.empty:
        return JSONResponse(status_code=400, content={"error": "The uploaded CSV contains no usable data."})

    try:
        columns = detect_columns(df)
        if taxon_column and taxon_column in df.columns:
            columns["taxon"] = taxon_column
        if sample_column and sample_column in df.columns:
            columns["sample"] = sample_column
        if site_column and site_column in df.columns:
            columns["site"] = site_column
        if abundance_column and abundance_column in df.columns:
            columns["abundance"] = abundance_column

        warnings = []
        quality = assess_quality(df, columns)
        abundance_mode = columns["abundance"] is not None
        if not abundance_mode:
            warnings.append("Abundance information was not available. Presence/absence analysis is being used.")

        taxonomy = analyze_taxonomy(df, columns)
        biodiversity = {"available": False, "reason": "No taxon column identified."}
        if columns["taxon"]:
            taxon_col = columns["taxon"]
            if abundance_mode:
                work = df[[taxon_col, columns["abundance"]]].copy()
                work[columns["abundance"]] = pd.to_numeric(work[columns["abundance"]], errors="coerce")
                work = work.dropna()
                counts = work.groupby(taxon_col)[columns["abundance"]].sum()
            else:
                taxa_series = df[taxon_col].dropna().astype(str).str.strip()
                taxa_series = taxa_series[taxa_series != ""]
                counts = taxa_series.value_counts()
            biodiversity = compute_biodiversity(counts, abundance_mode)

        sites = analyze_sites(df, columns, abundance_mode)
        community = compute_community_similarity(df, columns, abundance_mode)
        sequence_qc = analyze_sequence_qc(df, columns)
        spatial = analyze_spatial(df, columns)
        temporal = analyze_temporal(df, columns)

        result = {
            "is_demo": False,
            "dataset": {"rows": int(df.shape[0]), "columns": int(df.shape[1])},
            "columns": {
                "taxon": columns["taxon"], "sample": columns["sample"], "site": columns["site"],
                "abundance": columns["abundance"], "sequence": columns["sequence"],
                "latitude": columns["latitude"], "longitude": columns["longitude"], "date": columns["date"],
                "ranks": columns["ranks"], "all_columns": df.columns.tolist(),
                "confident": columns["taxon"] is not None,
                "analysis_mode": "abundance" if abundance_mode else "presence_absence",
            },
            "quality": quality,
            "taxonomy": taxonomy,
            "biodiversity": biodiversity,
            "sites": sites,
            "community": community,
            "sequence_qc": sequence_qc,
            "spatial": spatial,
            "temporal": temporal,
            "warnings": quality["warnings"] + warnings,
        }
        return JSONResponse(content=result)
    except Exception:
        return JSONResponse(
            status_code=500,
            content={"error": "An unexpected error occurred while analyzing the dataset with the provided column mapping."},
        )
