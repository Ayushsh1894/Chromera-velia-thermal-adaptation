# Chromera velia thermal-adaptation RNA-seq

Analysis code and processed data for Sharma *et al.*, *"Transcriptomic and
lipidomic profiles reveal greater thermal stability in a tropical versus a
temperate isolate of Chromera velia"*.

Two strains of *Chromera velia* — **SH** (Sydney Harbour, temperate, reference
strain CCMP2878) and **OTI** (One Tree Island, tropical, southern Great
Barrier Reef) — were grown at **20 °C** and **26 °C**, three biological
replicates per condition (12 RNA-seq libraries total). The analysis compares
the two strains at each temperature (strain difference), each strain's
response to temperature (thermal response), and the core/pan transcriptome.

> **A note on strain naming.** The manuscript and this repository use **OTI**
> throughout. A small number of legacy filenames inherited from the original
> CryptoDB exports may still read **OTIK** (e.g. `Sh_down_OTIK_up`) —
> these refer to the same strain.

---

## Repository layout

```
.
├── codes/                          R analysis scripts (run in numbered order)
├── excel sheets/                   annotated differential-expression + GO tables
│   ├── Strain Difference/
│   │   ├── OTI/                    genes higher in OTI (heat / cool)
│   │   └── SH/                     genes higher in SH  (heat / cool)
│   ├── Thermal response/
│   │   ├── OTI/                    OTI strain-unique responders (up / down)
│   │   ├── SH/                     SH strain-unique responders  (up / down)
│   │   └── concordant/             shared thermal DEGs by concordance category
│   └── experimental_design_temp.csv
└── climate data/                   in-situ sea-surface-temperature data for the
                                     two isolation sites (GBR vs Sydney)
```

> **Note:** the `Strain Difference` folder is currently misspelled
> (`Strain Diffrence`) in the repository — rename it to match the paths above.

---

## `codes/` — analysis scripts

R (≥ 4.2) with Bioconductor. Scripts are numbered in dependency order; run
from the folder containing the two input files below (or set `project_dir`
at the top of each script).

### Inputs

| File | Description |
|------|-------------|
| `Count_matrix.xlsx` | Gene × sample raw count matrix (30,492 genes × 12 samples), first column `Geneid`. Produced from raw reads by STAR alignment + featureCounts (see Methods). |
| `experimental_design_temp.csv` | Sample metadata: `sample`, `strain` (SH/OTI), `temperature` (T20/T26). Provided in `excel sheets/`. |

Raw sequence reads are deposited at NCBI under BioProject accession
**PRJNA1519432**. The count matrix is provided as supplementary data with the
article. The read-processing pipeline (Trimmomatic → STAR → featureCounts)
was run separately on an HPC cluster and is described in the Methods; this
repository begins from the count matrix.

### Scripts

| Script | What it does |
|--------|---------------|
| `01_core_transcriptome.R` | Core and pan transcriptome by presence/absence (edgeR CPM ≥ 1 in all 3 replicates); bar chart, 4-way Venn, UpSet. |
| `02_deseq2_master.R` | Builds the DESeq2 objects and all contrasts, saves `DEA_results.RData`; PCA and top-50 variable-gene heatmap. **Run this first** — the other scripts load its output. |
| `03_strain_comparison.R` | OTI-vs-SH at each temperature (per-temperature subsets); Venn; log2FC(26 °C)-vs-log2FC(20 °C) scatter. |
| `04_strain_directional_categories.R` | Splits strain DEGs into eight direction-resolved categories (constitutive / heat-specific / cool-specific / direction-flip, per strain); bar plot, UpSet, annotated quadrant scatter. Produces the strain-difference gene lists. |
| `05_thermal_concordance.R` | For genes DE in both strains' thermal responses, classifies concordant (same direction) vs discordant (opposite direction) genes; concordance scatter. |

### Key parameters

- **DE significance:** `padj < 0.05` **and** `|log2FoldChange| > 1` throughout.
- **Low-count filter:** keep genes with ≥ 10 counts in ≥ 3 samples.
- **Reference levels:** strain = SH, temperature = T20.
  - Strain contrasts: positive log2FC = higher in OTI.
  - Thermal contrasts: positive log2FC = up at 26 °C.

**Two normalisations, two purposes.** Presence/absence (script 01) uses
**edgeR CPM**; all differential expression (scripts 02–05) uses **DESeq2
median-of-ratios** normalisation. These answer different questions and are
not interchangeable.

**Within-strain thermal contrasts** are computed on **per-strain subsets**
(each `~ temperature`) rather than from an interaction term on the full
model, to avoid model-fitting issues arising from the very large strain
effect.

Each script installs any missing packages on first run and writes a
`sessionInfo_*.txt` recording exact package versions.

### Expected results (validation targets)

| Quantity | Value |
|----------|-------|
| Core transcriptome | 15,990 genes |
| Pan-transcriptome | 23,081 genes |
| Strain DEGs at 20 °C | 8,065 (5,284 OTI-higher / 2,781 SH-higher) |
| Strain DEGs at 26 °C | 8,253 (6,113 OTI-higher / 2,140 SH-higher) |
| Constitutive strain DEGs (both temperatures) | 6,455 |
| SH thermal DEGs (T26 vs T20) | 4,511 |
| OTI thermal DEGs (T26 vs T20) | 4,089 |
| Shared thermal DEGs | 1,614 |

---

## `excel sheets/` — annotated DE and GO tables

Gene lists exported by the scripts were submitted to
[CryptoDB](https://cryptodb.org) for functional annotation and Gene Ontology
enrichment. Each workbook contains three sheets with an identical structure:

| Sheet | Contents |
|-------|----------|
| `DEG` | The differentially expressed genes in that set, with log2 fold change and adjusted p-value. |
| `Annotation` | Gene IDs mapped to product descriptions (CryptoDB). |
| `Gene_ontology` | GO enrichment results (GO term, fold enrichment, p-value, Benjamini–Hochberg-corrected p-value), from CryptoDB's enrichment tool (Fisher's exact test). |

**Strain Difference** — genes differing between strains at a given
temperature (from `04_strain_directional_categories.R`):
- `OTI/Strain_OTI_heat_up`, `OTI/Strain_OTI_cool_up` — higher in OTI at 26 °C / 20 °C
- `SH/Strain_SH_heat_up`, `SH/Strain_SH_cool_up` — higher in SH at 26 °C / 20 °C

**Thermal response** — genes responding to temperature within a strain (from
the thermal contrasts and `05_thermal_concordance.R`):
- `OTI/OTI_unique_up`, `OTI/OTI_unique_down` — genes OTI uniquely up/down-regulates at 26 °C
- `SH/SH_unique_up`, `SH/SH_unique_down` — genes SH uniquely up/down-regulates at 26 °C
- `concordant/Concordant_up`, `concordant/Concordant_down` — up / down at 26 °C in both strains
- `concordant/Sh_down_OTIK_up`, `concordant/SH_up_OTI_down` — discordant (opposite direction between strains)

---

## `climate data/` — sea-surface temperature

In-situ sea-surface-temperature records (IMOS SOOP-SST) for the two isolation
sites, used to place the experimental temperatures (20 °C and 26 °C) in the
context of each strain's native thermal environment (temperate Sydney vs the
tropical Great Barrier Reef).

| File | Description |
|------|-------------|
| `01_raw_SOOP-SST_GBR.zip` | Raw SOOP-SST records, Great Barrier Reef |
| `02_raw_SOOP-SST_Sydney.csv` | Raw SOOP-SST records, Sydney |
| `03_cleaned_SOOP-SST_GBR.csv` | QC'd/cleaned GBR records |
| `04_cleaned_SOOP-SST_Sydney.csv` | QC'd/cleaned Sydney records |
| `05_climatology_GBR.csv` | Monthly climatology, GBR |
| `06_climatology_Sydney.csv` | Monthly climatology, Sydney |
| `07_annual_means_GBR.csv` | Annual mean SST, GBR |
| `08_annual_means_Sydney.csv` | Annual mean SST, Sydney |

Source and processing are described in the Methods.

---

## Environment

R (≥ 4.2) with Bioconductor. Main packages: **DESeq2**, **edgeR**,
**EnhancedVolcano**, **ggplot2**, **dplyr**, **tibble**, **ggvenn**,
**UpSetR**, **ggrepel**, **pheatmap**, **RColorBrewer**.

## Citation

If you use this code or data, please cite the article (details on
publication) and this repository.
