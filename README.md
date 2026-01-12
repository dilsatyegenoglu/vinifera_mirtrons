# vinifera_mirtrons

This repository contains the full implementation of the computational pipeline used in the manuscript submitted, it is an *in silico* discovery and classification of mirtron candidates in *Vitis vinifera*.
The pipeline is designed as a stepwise, modular framework reflecting the exact order, logic, and parameterization described in the **Materials and Methods** section of the manuscript.
## 1. Genome Annotation–Based Intron Extraction
Introns are extracted from exon annotations using genome FASTA and GFF files.
**Input**
- Genome FASTA (.fna)
- Genome annotation (.gff)
**Script**
- `1_intron_stat.py`
**Output**
- Raw intron sequences (FASTA)
- Intron coordinate table (TSV)
- Descriptive statistics (length distribution, strand, chromosome)
## 2. Intron Consolidation
Redundant introns originating from alternative transcript isoforms are collapsed into a non-redundant intron set based on genomic coordinates.
**Script**
- `2_consolidated_stat.py`
**Output**
- Consolidated intron FASTA/TSV
- Updated intron statistics
## 3. Length-Based Filtering
Introns are filtered based on length constraints to exclude excessively short or long introns that are incompatible with known mirtron characteristics.
**Script**
- `3_length_stat.py`
**Parameters**
- Minimum intron length (user-defined)
- Maximum intron length (user-defined)
**Output**
- Length-filtered intron set
- Length distribution statistics
## 4. GC Content Filtering
GC content is calculated for each intron and used as a filtering criterion.
**Script**
- `4_gcfilter_stat.py`
**Parameters**
- Minimum GC percentage
- Maximum GC percentage
**Output**
- GC-filtered intron set
- GC content distribution statistics
## 5. RNA Secondary Structure Prediction
Secondary structures are predicted using RNAfold to compute minimum free energy (MFE) values and dot-bracket structures.
**Script**
- `5_rnafold_stat.py`
**Requirements**
- ViennaRNA package (RNAfold)
**Output**
- Folded intron structures
- MFE, AMFE, and MFEI values
- Structural statistics
## 6. Structural Validation and Scoring
Each folded intron is evaluated for structural features associated with mirtrons, including stem length, bulge frequency, base-pairing, and overhang properties.
A composite confidence score is calculated and used to classify candidates into confidence levels.
**Script**
- `6_supervalidator_stat.py`
**Scoring Classes**
- `very_high`
- `high`
- `moderate`
- `low`
**Output**
- Classified mirtron candidates (FASTA/TSV)
- Structural feature tables
- Class-wise statistics
## 7. Sequence Similarity Analysis
Candidate mirtrons are compared against known precursor and/or mature miRNA datasets using BLAST.
**Script**
- `7_blast_stat.py`
**Output**
- BLAST match tables
- Known vs novel candidate annotation
## 8. Drosha-Free Processing Assessment
Candidates are evaluated for Drosha-independent processing based on 5′ and 3′ overhang constraints.
**Script**
- `8_overhang_stat.py`
**Criteria**
- Overhang length thresholds as defined in the manuscript
**Output**
- Drosha-free mirtron candidate subset
- Overhang statistics
## 9. Mirtron Category Subdivision
Mirtron candidates are subdivided into structural categories based on overhang orientation.
**Script**
- `9_mirtrondivider_stat.py`
**Categories**
- Canonical
- 5′-tailed
- 3′-tailed
- Both-tailed
**Output**
- Category-specific FASTA and TSV files
- Summary statistics
## Input Data Availability
Large genome and annotation files (e.g. `.fna`, `.gff`) are not included in this repository due to size constraints.  
They can be obtained from public databases such as NCBI Genome or Ensembl Plants.
## Reproducibility and Versioning
All pipeline steps are deterministic and parameter-driven.  
The exact version of the pipeline used in the manuscript is archived with a DOI via Zenodo and corresponds to a tagged GitHub release.
## License
This software is released under the MIT License.

Large genome and annotation files (e.g. `.fna`, `.gff`) are intentionally excluded from this repository and should be obtained from public databases.
## 📥 Input Data
Typical inputs include:
- Genome FASTA (`.fna`)
- Genome annotation (`.gff` or `.gtf`)
- Known plant mirtron datasets (precursor and/or mature FASTA)
Reference genome and annotation files for *Vitis vinifera* can be obtained from:
- NCBI Genome
- Ensembl Plants etc
## 📤 Output Data
The pipeline generates:
- Filtered intron FASTA and TSV files  
- Classified mirtron candidate lists  
- Category-specific FASTA/TSV outputs  
- Statistical summary tables and reports  

These outputs are designed to be directly usable for downstream analyses and manuscript preparation.
## 🔁 Reproducibility
All steps are deterministic and parameter-driven.  
Intermediate outputs are preserved to ensure transparency and reproducibility.
A tagged release corresponding to the version used in the associated manuscript is archived with a DOI.
## 📖 Citation

If you use this pipeline in your research, please cite the associated manuscript (in preparation) and reference the archived software version:

> vinifera_mirtrons: a structure-aware pipeline for mirtron discovery in *Vitis vinifera*.

(Full citation and DOI will be added upon publication.)

## 📬 Contact

For questions or collaboration inquiries, please contact:
  
**Dilsat Yegenoglu**  
Manisa Celal Bayar University, Alasehir Vocational School, Alasehir, Manisa, TURKIYE 
📧 dilsat.yegenoglu@cbu.edu.tr
** Tevfik Hasan Can**
Manisa Celal Bayar University, Alasehir Vocational School, Alasehir, Manisa, TURKIYE 
📧 tevfik.can@cbu.edu.tr

