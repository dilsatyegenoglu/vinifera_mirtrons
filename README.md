Mirtron Discovery & Validation Pipeline

A fully automated 9-module bioinformatics pipeline for genome-wide mirtron candidate discovery and validation.

The pipeline performs:

Intron extraction from genome + GFF

Length filtering

GC filtering

RNA secondary structure prediction (RNAfold)

Structural validation

Data-driven BLAST optimization (autoblast3.py)

Drosha-free overhang classification

Final mirtron categorization

📋 Table of Contents

Requirements

Installation

Input Data

Workflow Overview

Module Descriptions

BLAST Optimization (autoblast3.py)

Output Structure

Troubleshooting

Citation

🔧 Requirements
System

Python ≥ 3.8

RAM ≥ 8 GB (16 GB recommended)

Linux / WSL2 recommended

External Tools
Tool	Test Command	Purpose
RNAfold (ViennaRNA)	RNAfold --version	Secondary structure prediction
BLAST+	blastn -version	Similarity search
makeblastdb	makeblastdb -version	BLAST database creation

Install (Ubuntu/WSL):

sudo apt update
sudo apt install viennarna ncbi-blast+

Python Packages
pip install biopython numpy pandas matplotlib


Plotly is optional (for interactive plots):

pip install plotly kaleido

📦 Installation
git clone https://github.com/USERNAME/mirtron-discovery-pipeline.git
cd mirtron-discovery-pipeline

📁 Required Input Files
File	Format	Description
genome.fa	FASTA	Reference genome
annotation.gff	GFF3	Genome annotation
premature.fasta	FASTA	miRNA precursor database
mature.fasta	FASTA	miRNA mature database

⚠ FASTA and GFF chromosome names must match exactly.

🚀 Workflow Overview
Step	Script	Description
1	1_intron_stat.py	Intron extraction
2	2_consolidated_stat.py	Consolidation
3	3_length_stat.py	Length filtering
4	4_gcfilter_stat.py	GC filtering
5	5_rnafold_stat.py	RNAfold analysis
6	6_supervalidator_stat.py	Structural validation
7	autoblast3.py	Data-driven BLAST optimization
8	8_overhang_stat.py	Drosha-free classification
9	9_mirtrondivider_stat.py	Final categorization
📌 Module Descriptions
🟢 Module 1 – Intron Extraction
python3 1_intron_stat.py


Output:

intron_outputs/

🟢 Module 2 – Consolidation
python3 2_consolidated_stat.py


Output:

consolidated_introns/

🟢 Module 3 – Length Filtering
python3 3_length_stat.py


User input:

Minimum length

Maximum length

Output:

length_filter_results/

🟢 Module 4 – GC Filtering
python3 4_gcfilter_stat.py


User input:

Minimum GC%

Maximum GC%

Output:

gc_filter_results/

🟢 Module 5 – RNAfold
python3 5_rnafold_stat.py


User input:

Minimum MFE

Maximum MFE

Output:

rnafold_results/

🟢 Module 6 – Structural Validator
python3 6_supervalidator_stat.py


Creates:

validator_results/
    └── classes/
        ├── very_high.fa
        ├── high.fa
        ├── moderate.fa
        └── low.fa

🚀 Module 7 – BLAST Optimization (autoblast3.py)

Script:

python3 autoblast3.py


This is NOT a simple BLAST wrapper.

It performs:

Exploratory BLAST (broad search)

Data-driven threshold optimization

Percentile-based identity selection

E-value MIN & MAX optimization

Identity–Coverage correlation analysis

Combined score calculation

Cross-validation

Category-specific thresholding

Maximum 5 hits per intron (highest bitscore)

Automatic report generation

Three modes:

Full automatic (recommended)

Optimization only

Manual thresholds

Output:

blast_results_ultimate/
    ├── premature/
    ├── mature/
    ├── optimization/
    ├── validation/
    └── summary.txt


This replaces the older 7_blast_cov_1.py.

🟢 Module 8 – Drosha-Free Overhang
python3 8_overhang_stat.py


User selects maximum allowed overhang (1–3 nt).

Output:

drosha_free_results/

🟢 Module 9 – Final Categorization
python3 9_mirtrondivider_stat.py


Final output:

mirtron_divide/

📂 Final Output Structure
mirtron_divide/                ← Final categorized mirtrons
blast_results_ultimate/        ← Optimized BLAST results
drosha_free_results/           ← Overhang classification
validator_results/             ← Structural classes
rnafold_results/               ← MFE-filtered sequences

❓ Troubleshooting

RNAfold not found:

sudo apt install viennarna


BLAST not found:

sudo apt install ncbi-blast+


Missing Python module:

pip install <module_name>


👩‍🔬 Developers

Dr. Emine Dilşat Yegenoğlu
Hasan Tevfik Can

Manisa Celal Bayar University
DNA & RNA Regulators Bioinformatics Laboratory

Version: 2.0
Last Updated: February 2026