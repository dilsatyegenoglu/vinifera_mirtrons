#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Mirtron Discovery Script – Step-by-step
# Authors: Emine Dilsat Yegenoglu & Tevfik Hasan Can

print("""REGAVCI: MODULE I, Hunt for the Mirtrons""")

print("""
Mirtron Analysis Script
========================
This script performs RNA secondary structure prediction using RNAfold
on intron sequences and classifies them based on dot-bracket structure and MFE values.
It is part of a larger pipeline for mirtron discovery and filtering.

⚠️ WARNING LIST:
---------------
1. GFF and genome FASTA files must contain matching chromosome names.
   Example: 'Chr1' is NOT equal to 'chr1'.

2. All FASTA files must be properly formatted:
   - Each sequence must start with a '>' symbol.
   - No blank lines should exist between sequences.

3. RNAfold must be installed and available in the system PATH.
   - Test with: `RNAfold --version`

4. Intron sequences should be at least ~30 nt long.

5. Avoid using non-ASCII characters in file names.

6. Extremely long sequences may timeout during RNAfold execution.

7. Dot-bracket results with no clear hairpin may be classified as 'ambiguous'.

Dependencies:
-------------
- RNAfold (ViennaRNA)
- Python ≥ 3.6
- Biopython
""")

print("""
# Mirtron Discovery Script – Step-by-step
# Developers: Assist. Prof. Dr. Emine Dilsat Yegenoglu & Lecturer Tevfik Hasan Can, M.Sc.
# Manisa Celal Bayar University, Alasehir Vocational School, Department of Plant and Animal Sciences, DNA & RNA Regulators Bioinformatics Laboratory.
# Collaborators: 
# Prof. Dr. Murat Afsar, Manisa Celal Bayar University, Faculty of Engineering and Natural Sciences, Department of Biology, Unit of Zoology
# Prof. Dr. Meltem Sesli, Manisa Celal Bayar University, School of Tobacco Expertise, Department of Tobacco Quality and Breeding
# Prof. Dr. Yakut Gevrekci, Ege University, Agriculture Faculty, Department of Animal Science, Unit of Biometry-Genetics
# Assist. Prof. Dr. Emine Berberoglu, Tokat Gaziosmanpasa University, Agriculture Faculty, Department of Animal Science, Unit of Biometry-Genetics
# Assist. Prof. Dr. Tamer Kusaksiz, Manisa Celal Bayar University, Alasehir Vocational School, Department of Plant and Animal Sciences
""")
print("""

STEP 1: INTRON EXTRACTION
==========================
Extracts introns from GFF + Genome FASTA
Output: intron_outputs/introns.tsv & introns.fa
Enhanced with comprehensive statistics
""")

import sys, os
from collections import defaultdict, namedtuple, Counter
import statistics
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple

try:
    from pyfaidx import Fasta
except ImportError:
    print("[ERROR] pyfaidx is not installed. Run: pip install pyfaidx")
    sys.exit(1)

Feature = namedtuple("Feature", "seqid source ftype start end score strand phase attrs raw")

def ask(prompt: str, default: str = None) -> str:
    try:
        s = input(prompt)
        if not s and default is not None:
            return default
        return s
    except KeyboardInterrupt:
        print("\n[ABORT] User interrupted.")
        sys.exit(1)

def parse_attributes(attr: str) -> Dict[str, str]:
    d = {}
    for item in attr.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def parse_gff(path: str, allowed_seqids: set) -> List[Feature]:
    feats = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            if allowed_seqids and seqid not in allowed_seqids:
                continue
            try:
                s = int(start); e = int(end)
            except:
                continue
            feats.append(Feature(seqid, source, ftype, s, e, score, strand, phase, parse_attributes(attrs), line))
    return feats

_RC_TABLE = str.maketrans({
    'A':'T','C':'G','G':'C','T':'A','U':'A',
    'a':'t','c':'g','g':'c','t':'a','u':'a',
    'N':'N','n':'n'
})
def revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]

def load_fasta_with_headers(fasta_path: str) -> Tuple[Fasta, Dict[str, Tuple[int,str]]]:
    fa = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    info = {}
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(">"):
                head = line[1:].strip()
                parts = head.split(None, 1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                info[name] = (0, desc)
    for name in list(info.keys()):
        try:
            info[name] = (len(fa[name]), info[name][1])
        except KeyError:
            pass
    return fa, info

def fetch_seq(fa: Fasta, seqid: str, start: int, end: int, strand: str) -> str:
    s = min(start, end)
    e = max(start, end)
    seq = str(fa[seqid][s-1:e])
    if strand == "-":
        seq = revcomp(seq)
    return seq

def introns_from_exons(exon_list: List[Feature]) -> List[Tuple[int,int]]:
    if len(exon_list) < 2:
        return []
    ex_sorted = sorted(exon_list, key=lambda x: x.start)
    intrs = []
    for i in range(len(ex_sorted)-1):
        left = ex_sorted[i]
        right = ex_sorted[i+1]
        intr_start = left.end + 1
        intr_end   = right.start - 1
        if intr_start <= intr_end:
            intrs.append((intr_start, intr_end))
    return intrs

def print_contig_table(info: Dict[str, Tuple[int,str]]) -> None:
    print("\n📜 Available contigs/chromosomes:")
    print("No\tName\tLength\tDescription")
    for idx, name in enumerate(info.keys(), start=1):
        length, desc = info[name]
        print(f"{idx}\t{name}\t{length}\t{desc}")

def calculate_statistics(intron_lengths, introns_per_transcript, strand_counts, seqid_counts):
    """Calculate comprehensive statistics for introns"""
    stats = {}
    
    # Basic length statistics
    stats['total_introns'] = len(intron_lengths)
    stats['mean_length'] = statistics.mean(intron_lengths) if intron_lengths else 0
    stats['median_length'] = statistics.median(intron_lengths) if intron_lengths else 0
    stats['min_length'] = min(intron_lengths) if intron_lengths else 0
    stats['max_length'] = max(intron_lengths) if intron_lengths else 0
    stats['std_length'] = statistics.stdev(intron_lengths) if len(intron_lengths) > 1 else 0
    
    # Percentiles
    if intron_lengths:
        stats['percentiles'] = {
            p: np.percentile(intron_lengths, p) 
            for p in [25, 50, 75, 90, 95, 99]
        }
    else:
        stats['percentiles'] = {}
    
    # Introns per transcript
    stats['mean_introns_per_transcript'] = statistics.mean(introns_per_transcript) if introns_per_transcript else 0
    stats['max_introns_per_transcript'] = max(introns_per_transcript) if introns_per_transcript else 0
    stats['transcript_count'] = len(introns_per_transcript)
    
    # Strand distribution
    stats['strand_plus'] = strand_counts.get('+', 0)
    stats['strand_minus'] = strand_counts.get('-', 0)
    stats['strand_unknown'] = strand_counts.get('.', 0) + strand_counts.get('?', 0)
    
    # Chromosome/contig distribution
    stats['seqid_distribution'] = dict(seqid_counts.most_common())
    stats['unique_seqids'] = len(seqid_counts)
    
    # Length categories
    stats['length_categories'] = {
        'very_short': sum(1 for l in intron_lengths if l < 50),
        'short': sum(1 for l in intron_lengths if 50 <= l < 100),
        'medium': sum(1 for l in intron_lengths if 100 <= l < 500),
        'long': sum(1 for l in intron_lengths if 500 <= l < 2000),
        'very_long': sum(1 for l in intron_lengths if l >= 2000)
    }
    
    return stats

def generate_plots(intron_lengths, strand_counts, seqid_counts, outdir):
    """Generate comprehensive statistical plots"""
    plot_dir = os.path.join(outdir, "statistics_plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Length distribution histogram
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.hist(intron_lengths, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('Intron Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Intron Length Distribution')
    plt.yscale('log')
    
    # 2. Strand distribution
    plt.subplot(2, 2, 2)
    strands = ['+', '-', 'Unknown']
    counts = [strand_counts.get('+', 0), strand_counts.get('-', 0), 
              strand_counts.get('.', 0) + strand_counts.get('?', 0)]
    colors = ['lightgreen', 'lightcoral', 'lightgray']
    plt.pie(counts, labels=strands, autopct='%1.1f%%', colors=colors)
    plt.title('Strand Distribution')
    
    # 3. Top 10 chromosomes/contigs
    plt.subplot(2, 2, 3)
    top_seqids = seqid_counts.most_common(10)
    if top_seqids:
        names, counts = zip(*top_seqids)
        plt.bar(range(len(names)), counts, color='orange', alpha=0.7)
        plt.xticks(range(len(names)), names, rotation=45, ha='right')
        plt.xlabel('Chromosome/Contig')
        plt.ylabel('Intron Count')
        plt.title('Top 10 Chromosomes/Contigs')
    
    # 4. Length categories
    plt.subplot(2, 2, 4)
    categories = ['<50', '50-100', '100-500', '500-2000', '≥2000']
    category_counts = [
        sum(1 for l in intron_lengths if l < 50),
        sum(1 for l in intron_lengths if 50 <= l < 100),
        sum(1 for l in intron_lengths if 100 <= l < 500),
        sum(1 for l in intron_lengths if 500 <= l < 2000),
        sum(1 for l in intron_lengths if l >= 2000)
    ]
    plt.bar(categories, category_counts, color='purple', alpha=0.7)
    plt.xlabel('Length Category (bp)')
    plt.ylabel('Count')
    plt.title('Intron Length Categories')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'intron_statistics_overview.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Detailed length histogram (log scale)
    plt.figure(figsize=(10, 6))
    plt.hist(intron_lengths, bins=100, alpha=0.7, color='steelblue', edgecolor='black')
    plt.xlabel('Intron Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Detailed Intron Length Distribution (Log Scale)')
    plt.yscale('log')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'length_distribution_detailed.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[STATS] Plots saved to {plot_dir}/")

def save_statistics_report(stats, outdir):
    """Save comprehensive statistics report"""
    report_path = os.path.join(outdir, "extraction_statistics_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("INTRON EXTRACTION STATISTICS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("SUMMARY STATISTICS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total introns extracted: {stats['total_introns']:,}\n")
        f.write(f"Total transcripts: {stats['transcript_count']:,}\n")
        f.write(f"Mean introns per transcript: {stats['mean_introns_per_transcript']:.2f}\n")
        f.write(f"Max introns per transcript: {stats['max_introns_per_transcript']}\n\n")
        
        f.write("LENGTH STATISTICS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Mean length: {stats['mean_length']:.2f} bp\n")
        f.write(f"Median length: {stats['median_length']:.2f} bp\n")
        f.write(f"Min length: {stats['min_length']} bp\n")
        f.write(f"Max length: {stats['max_length']} bp\n")
        f.write(f"Standard deviation: {stats['std_length']:.2f} bp\n\n")
        
        f.write("LENGTH PERCENTILES:\n")
        f.write("-" * 40 + "\n")
        for p, value in stats['percentiles'].items():
            f.write(f"{p}th percentile: {value:.2f} bp\n")
        f.write("\n")
        
        f.write("LENGTH CATEGORIES:\n")
        f.write("-" * 40 + "\n")
        cats = stats['length_categories']
        f.write(f"Very short (<50 bp): {cats['very_short']:,} ({cats['very_short']/stats['total_introns']*100:.2f}%)\n")
        f.write(f"Short (50-100 bp): {cats['short']:,} ({cats['short']/stats['total_introns']*100:.2f}%)\n")
        f.write(f"Medium (100-500 bp): {cats['medium']:,} ({cats['medium']/stats['total_introns']*100:.2f}%)\n")
        f.write(f"Long (500-2000 bp): {cats['long']:,} ({cats['long']/stats['total_introns']*100:.2f}%)\n")
        f.write(f"Very long (≥2000 bp): {cats['very_long']:,} ({cats['very_long']/stats['total_introns']*100:.2f}%)\n\n")
        
        f.write("STRAND DISTRIBUTION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Plus strand (+): {stats['strand_plus']:,} ({stats['strand_plus']/stats['total_introns']*100:.2f}%)\n")
        f.write(f"Minus strand (-): {stats['strand_minus']:,} ({stats['strand_minus']/stats['total_introns']*100:.2f}%)\n")
        f.write(f"Unknown strand: {stats['strand_unknown']:,} ({stats['strand_unknown']/stats['total_introns']*100:.2f}%)\n\n")
        
        f.write("CHROMOSOME/CONTIG DISTRIBUTION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Unique chromosomes/contigs: {stats['unique_seqids']}\n")
        for seqid, count in list(stats['seqid_distribution'].items())[:15]:  # Top 15
            f.write(f"  {seqid}: {count:,} ({count/stats['total_introns']*100:.2f}%)\n")
        if len(stats['seqid_distribution']) > 15:
            f.write(f"  ... and {len(stats['seqid_distribution']) - 15} more\n")
    
    print(f"[STATS] Statistics report saved to {report_path}")

def main():
    print("=" * 60)
    print("STEP 1: INTRON EXTRACTION (ENHANCED STATISTICS)")
    print("=" * 60)
    
    fasta_path = ask("Genome FASTA filename (path): ").strip()
    gff_path   = ask("GFF3 filename (path): ").strip()

    if not os.path.exists(fasta_path):
        print(f"[ERROR] FASTA not found: {fasta_path}")
        sys.exit(1)
    if not os.path.exists(gff_path):
        print(f"[ERROR] GFF not found: {gff_path}")
        sys.exit(1)

    fa, fa_info = load_fasta_with_headers(fasta_path)
    print_contig_table(fa_info)

    ans = ask("\nDo you want to exclude any contigs? (y/n) [n]: ", default="n").strip().lower()
    excluded_idx = set()
    excluded_names = []
    if ans in ("y", "yes"):
        raw = ask("Enter indices to exclude (comma-separated, e.g., 2,5,7): ").strip()
        if raw:
            for tok in raw.split(","):
                tok = tok.strip()
                if tok.isdigit():
                    excluded_idx.add(int(tok))
    ordered_names = list(fa_info.keys())
    allowed = set()
    for i, name in enumerate(ordered_names, start=1):
        if i not in excluded_idx:
            allowed.add(name)
        else:
            excluded_names.append(name)

    if not allowed:
        print("[ERROR] All contigs excluded. Exiting.")
        sys.exit(1)

    print("\n[Reading and filtering GFF ...]")
    feats = parse_gff(gff_path, allowed_seqids=allowed)
    if not feats:
        print("[ERROR] No features matched the selected contigs.")
        sys.exit(1)

    transcripts = defaultdict(dict)
    exons_all = []

    def get_parent(f):
        a = f.attrs
        return a.get("Parent") or a.get("parent") or a.get("transcript_id")

    print("[Indexing features ...]")
    for f in feats:
        t = f.ftype.lower()
        if t == "exon":
            exons_all.append(f)
            pid = get_parent(f)
            if not pid: continue
            if "seqid" not in transcripts[pid]:
                transcripts[pid]["seqid"] = f.seqid
            transcripts[pid]["strand"] = f.strand
            transcripts[pid].setdefault("exons", []).append(f)

    outdir = "intron_outputs"
    os.makedirs(outdir, exist_ok=True)

    intron_fasta = open(os.path.join(outdir, "introns.fa"), "w")
    intron_tsv   = open(os.path.join(outdir, "introns.tsv"), "w", encoding="utf-8")
    intron_tsv.write("ID\tseqid\tstart\tend\tstrand\tlength\tparent_transcript\tsequence\n")

    # Statistics collection
    intron_lengths = []
    introns_per_transcript = []
    strand_counts = Counter()
    seqid_counts = Counter()
    
    intron_count = 0
    print("\n[Intron extraction begins]")
    for tid, rec in transcripts.items():
        exs = rec.get("exons", [])
        if not exs or len(exs) < 2:
            continue
        seqid = rec.get("seqid")
        strand = rec.get("strand", "+")
        intrs = introns_from_exons(exs)
        
        # Record introns per transcript
        introns_per_transcript.append(len(intrs))
        
        for idx, (s, e) in enumerate(intrs, start=1):
            seq = fetch_seq(fa, seqid, s, e, strand)
            length = len(seq)
            iid = f"{seqid}:{s}-{e}:{strand}|{tid}|intron{idx}"
            intron_fasta.write(f">{iid}\n{seq}\n")
            intron_tsv.write(f"{iid}\t{seqid}\t{s}\t{e}\t{strand}\t{length}\t{tid}\t{seq}\n")
            
            # Collect statistics
            intron_lengths.append(length)
            strand_counts[strand] += 1
            seqid_counts[seqid] += 1
            
            intron_count += 1
            if intron_count % 1000 == 0:
                print(f"[intron] {intron_count:,} introns extracted...", flush=True)
    
    intron_fasta.close()
    intron_tsv.close()

    print(f"\n[EXTRACTION COMPLETED] Total introns: {intron_count:,}")
    
    # Generate statistics and plots
    if intron_count > 0:
        print("\n[STATS] Generating comprehensive statistics...")
        stats = calculate_statistics(intron_lengths, introns_per_transcript, strand_counts, seqid_counts)
        save_statistics_report(stats, outdir)
        generate_plots(intron_lengths, strand_counts, seqid_counts, outdir)
        
        # Terminal summary
        print("\n" + "=" * 60)
        print("EXTRACTION SUMMARY")
        print("=" * 60)
        print(f"Total introns: {stats['total_introns']:,}")
        print(f"Mean length: {stats['mean_length']:.1f} bp")
        print(f"Median length: {stats['median_length']:.1f} bp")
        print(f"Length range: {stats['min_length']:,} - {stats['max_length']:,} bp")
        print(f"Strand distribution: +:{stats['strand_plus']:,} (-:{stats['strand_minus']:,})")
        print(f"Chromosomes/contigs: {stats['unique_seqids']}")
        print(f"Transcripts processed: {stats['transcript_count']:,}")
        print(f"Mean introns/transcript: {stats['mean_introns_per_transcript']:.2f}")
        print("=" * 60)
    
    print(f"\n[STEP 1 COMPLETED]")
    print(f"Output files: {outdir}/introns.tsv, {outdir}/introns.fa")
    print(f"Statistics: {outdir}/extraction_statistics_report.txt")
    print(f"Plots: {outdir}/statistics_plots/")

if __name__ == "__main__":
    main()