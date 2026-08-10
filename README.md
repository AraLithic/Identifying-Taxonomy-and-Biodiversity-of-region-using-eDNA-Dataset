# TerraTrace by Sayed Umair Ali

**Deterministic eDNA & Biodiversity Analysis**

TerraTrace is a local web application for analyzing processed **eDNA / metabarcoding datasets** built for SIH. It provides deterministic taxonomy, biodiversity, community, sequence-quality, spatial, and temporal analyses directly from uploaded CSV datasets.

The application runs entirely locally with a lightweight **FastAPI backend** and a static **HTML/JavaScript frontend**.

> **No AI · No database · No external analysis APIs used**

---

## Features

* CSV-based eDNA / metabarcoding dataset analysis
* Automatic dataset column detection
* Manual column mapping when automatic detection is uncertain
* Taxonomic classification and summaries
* Biodiversity metrics
* Site-level comparisons
* Community similarity analysis
* Sequence quality control
* Spatial analysis
* Temporal analysis
* Synthetic demonstration dataset
* Interactive visualizations with Chart.js
* Local processing with no database or external analysis service

TerraTrace processes uploaded datasets in memory using Python scientific-computing libraries.

---

## Project Structure

```text
TerraTrace/
├── backend/
│   └── app.py
│
├── frontend/
│   └── index.html
│
└── README.md
```

---

## Requirements

* **Python 3.9+**
* **Git**
* A modern web browser

### Python dependencies

```text
fastapi
uvicorn
pandas
numpy
scipy
python-multipart
```

---

## Installation

### 1. Clone the repository

```bash
git clone <REPOSITORY_URL>
cd TerraTrace
```

Replace `<REPOSITORY_URL>` with the URL of this repository.

---

### 2. Create a virtual environment

```bash
python3 -m venv venv
```

Activate it:

#### macOS / Linux

```bash
source venv/bin/activate
```

#### Windows

```powershell
venv\Scripts\activate
```

---

### 3. Install dependencies

```bash
pip install fastapi uvicorn pandas numpy scipy python-multipart
```

---

## Running Locally

TerraTrace uses two local servers:

* **Backend:** `http://127.0.0.1:8000`
* **Frontend:** `http://127.0.0.1:5500`

Both servers must be running at the same time.

### Terminal 1 — Backend

From the project root:

```bash
source venv/bin/activate
uvicorn backend.app:app --reload
```

The FastAPI server will start on:

```text
http://127.0.0.1:8000
```

The backend exposes health, demo, analysis, and manual-mapping endpoints.

### Terminal 2 — Frontend

Open a second terminal:

```bash
cd TerraTrace
python3 -m http.server 5500 --directory frontend
```

Then open:

**http://127.0.0.1:5500**

The frontend is configured to communicate with the backend at `127.0.0.1:8000`.

---

## Quick Start

Once the repository has been cloned:

### Terminal 1

```bash
cd TerraTrace
python3 -m venv venv
source venv/bin/activate
pip install fastapi uvicorn pandas numpy scipy python-multipart
uvicorn backend.app:app --reload
```

### Terminal 2

```bash
cd TerraTrace
python3 -m http.server 5500 --directory frontend
```

Open:

```text
http://127.0.0.1:5500
```

---

## Using TerraTrace

1. Open the frontend in your browser.
2. Upload a `.csv` eDNA / metabarcoding dataset.
3. Click **Analyze Dataset**.
4. Review the generated analysis.
5. If column detection is ambiguous, use the **Column Mapping** interface and re-run the analysis.

The frontend accepts CSV files and provides a synthetic demo dataset for testing the application without an input file.

---

## Input Data

TerraTrace automatically detects common column names for key dataset fields.

| Field     | Supported column names                            |
| --------- | ------------------------------------------------- |
| Taxon     | `species`, `scientific_name`, `taxon`, `organism` |
| Sample    | `sample_id`, `sample`, `sample_name`              |
| Site      | `site`, `location`, `sampling_site`, `station`    |
| Abundance | `reads`, `read_count`, `abundance`, `count`       |
| Sequence  | `sequence`, `seq`, `dna_sequence`                 |
| Latitude  | `latitude`, `lat`                                 |
| Longitude | `longitude`, `lon`, `lng`, `long`                 |
| Date      | `date`, `sampling_date`, `collection_date`        |

Taxonomic rank columns such as `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, and `species` are also recognized.

If a suitable taxon column cannot be confidently identified, TerraTrace allows the user to manually map the relevant columns before re-analysis.

---

## Analysis

TerraTrace provides analysis across several areas:

### Dataset Quality

* Row and column counts
* Missing values
* Duplicate rows
* Invalid abundance values
* Invalid coordinates
* Invalid sequences
* Dataset warnings

### Taxonomy

* Taxonomic summaries
* Kingdom through species-level information
* Observed taxa

### Biodiversity

* Diversity metrics
* Richness metrics
* Evenness
* Site-level biodiversity comparisons

### Community Analysis

* Site comparisons
* Community similarity
* Presence/absence analysis
* Abundance-based analysis when abundance data is available

### Sequence Quality Control

* Sequence validation
* Sequence length statistics
* GC-related analysis

### Spatial Analysis

* Geographic coordinates
* Site-level spatial information
* Spatial comparisons where sufficient data is available

### Temporal Analysis

* Date-based observations
* Temporal summaries and comparisons

The frontend organizes these results into dedicated analysis sections, including taxonomy, biodiversity, sites, community similarity, sequence QC, spatial analysis, and temporal analysis.

---

## API

The backend is implemented using FastAPI.

### Health Check

```http
GET /api/health
```

Returns the current backend status.

### Demo Dataset

```http
GET /api/demo
```

Runs TerraTrace against its built-in synthetic dataset.

### Analyze Dataset

```http
POST /api/analyze
```

Accepts a CSV file for analysis.

### Analyze With Manual Mapping

```http
POST /api/analyze-mapped
```

Accepts a CSV file together with optional manual mappings for:

```text
taxon_column
sample_column
site_column
abundance_column
```

---

## Upload Limit

The current maximum CSV upload size is **25 MB**.

This limit is defined in the FastAPI backend.

---

## Architecture

```text
                 ┌────────────────────────┐
                 │        Browser         │
                 │                        │
                 │   TerraTrace Frontend  │
                 │    frontend/index.html │
                 └───────────┬────────────┘
                             │
                         HTTP :5500
                             │
                             ▼
                 ┌────────────────────────┐
                 │   Python HTTP Server   │
                 │      Static Files      │
                 └────────────────────────┘


                 Browser API Requests
                             │
                         HTTP :8000
                             │
                             ▼
                 ┌────────────────────────┐
                 │        FastAPI         │
                 │     backend/app.py     │
                 └───────────┬────────────┘
                             │
                             ▼
                 ┌────────────────────────┐
                 │  pandas / NumPy /      │
                 │        SciPy           │
                 │                        │
                 │   CSV → Analysis       │
                 └────────────────────────┘
```

All dataset analysis is performed locally by the backend. The project does not require a database or external analysis API.

---

## Development

For backend development, use:

```bash
uvicorn backend.app:app --reload
```

The `--reload` option automatically restarts the backend when Python source files are changed.

For frontend development, keep the static server running:

```bash
python3 -m http.server 5500 --directory frontend
```

After modifying `frontend/index.html`, refresh the browser.

---

## Troubleshooting

### Backend is unreachable

Verify that the backend is running:

```bash
uvicorn backend.app:app --reload
```

Then visit:

```text
http://127.0.0.1:8000/api/health
```

Expected response:

```json
{
  "status": "ok",
  "service": "TerraTrace"
}
```

---

### Python module not found

Activate the virtual environment:

```bash
source venv/bin/activate
```

Then reinstall the dependencies:

```bash
pip install fastapi uvicorn pandas numpy scipy python-multipart
```

---

### Frontend does not load

Ensure the frontend server is running:

```bash
python3 -m http.server 5500 --directory frontend
```

Then open:

```text
http://127.0.0.1:5500
```

---

### Charts are unavailable

## TerraTrace loads Chart.js from a CDN and includes a fallback source. If Chart.js cannot be loaded, numerical analysis can still be displayed while charts are skipped.

## Stopping the Application

Press:

```text
Ctrl + C
```

in both terminals.

To deactivate the Python virtual environment:

```bash
deactivate
```

---

## Data & Privacy

TerraTrace is designed for local processing.

Uploaded CSV data is sent from the browser to the locally running FastAPI backend and processed in memory. The backend does not use a database or external analysis API.

---

## License

```text
MIT License
```

---

## Acknowledgements

Built for deterministic analysis of processed **eDNA / metabarcoding datasets** with a focus on taxonomy, biodiversity, community structure, sequence quality, spatial patterns, and temporal patterns.
