# Repository Explanation and Reproducibility Notes

This document provides a detailed clarification of the pipeline
structure, execution logic, parameter handling, and reproducibility
design of the Mirtron Discovery & Validation Pipeline.

------------------------------------------------------------------------

# 1. Full Pipeline Execution (Step-by-Step)

After placing the required input files in the project root directory:

-   genome.fa
-   annotation.gff
-   premature.fasta
-   mature.fasta

Run the modules sequentially:

``` bash
python3 1_intron_stat.py
python3 2_consolidated_stat.py
python3 3_length_stat.py
python3 4_gcfilter_stat.py
python3 5_rnafold_stat.py
python3 6_supervalidator_stat.py
python3 autoblast3.py
python3 8_overhang_stat.py
python3 9_mirtrondivider_stat.py
```

Final categorized mirtron candidates will be located in:

    mirtron_divide/

Each module produces an output directory that serves as input for the
next module.

------------------------------------------------------------------------

# 2. Parameter Handling

## Length Filtering

User-defined minimum and maximum nucleotide length.

## GC Filtering

User-defined GC percentage range.

## RNAfold (Structural Filtering)

User-defined Minimum Free Energy (MFE) range.

## BLAST Optimization (autoblast3.py)

The BLAST module performs:

-   Exploratory BLAST (broad threshold)
-   Percentile-based identity estimation
-   Empirical E-value range determination
-   Identity--coverage correlation analysis
-   Combined score calculation
-   Cross-validation (if sufficient data)
-   Category-specific threshold refinement
-   Maximum 5 hits per intron (highest bitscore)

All thresholds are derived empirically from the observed dataset
distribution.

------------------------------------------------------------------------

# 3. Repository Structure

    mirtron-discovery-pipeline/
    │
    ├── 1_intron_stat.py
    ├── 2_consolidated_stat.py
    ├── 3_length_stat.py
    ├── 4_gcfilter_stat.py
    ├── 5_rnafold_stat.py
    ├── 6_supervalidator_stat.py
    ├── autoblast3.py
    ├── 8_overhang_stat.py
    ├── 9_mirtrondivider_stat.py
    │
    ├── README.md
    ├── requirements.txt
    └── example_data/ (optional for testing)

------------------------------------------------------------------------

# 4. Data Flow Logic

1.  Genome + GFF → intron extraction\
2.  Intron FASTA → consolidation\
3.  Consolidated introns → length filtering\
4.  Length-filtered introns → GC filtering\
5.  GC-filtered introns → RNA secondary structure analysis\
6.  RNAfold-filtered introns → structural validation\
7.  Validated candidates → BLAST similarity analysis\
8.  BLAST-supported candidates → overhang classification\
9.  Overhang-filtered candidates → final mirtron categorization

Each step preserves traceability through TSV and FASTA outputs.

------------------------------------------------------------------------

# 5. Determinism and Reproducibility

-   No stochastic algorithms are used.
-   No random sampling is applied during filtering.
-   Thresholds are computed deterministically from input data.
-   Running the pipeline with identical input files produces identical
    outputs.
-   BLAST hit limitation (maximum 5 hits per intron) is deterministic
    and based on bitscore ranking.

------------------------------------------------------------------------

# 6. Intended Use

This repository is designed for:

-   Genome-wide mirtron candidate discovery
-   Structural validation of intronic hairpins
-   Empirical BLAST-supported similarity assessment
-   Reproducible computational regulatory RNA research

------------------------------------------------------------------------

# 7. Versioning

Version: 2.1\
Last Updated: February 2026

------------------------------------------------------------------------

End of document.
