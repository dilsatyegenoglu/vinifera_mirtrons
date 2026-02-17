#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RNAFOLD_ENHANCED SMART PARSER + STATISTICS
=======================================================
- Safe FASTA header parser (handles ANY number of | and =)
- Never crashes on malformed headers
- Pipeline-safe CoordID extraction
- Realtime RNAfold with ETA
- MFE range filtering
- IN/OUT FASTA + TSV generation
- Enhanced statistical analysis
- Comprehensive visualizations
- Terminal + file summary
- All original outputs preserved
"""

import os
import subprocess
import time
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import statistics
import numpy as np
from collections import Counter

# ================================================================
# PATHS
# ================================================================
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
GC_DIR   = os.path.join(BASE_DIR, "gc_filter_results")
IN_FA    = os.path.join(GC_DIR, "gc_inrange.fa")

OUT_DIR  = os.path.join(BASE_DIR, "rnafold_results")
os.makedirs(OUT_DIR, exist_ok=True)

PLOT_DIR = os.path.join(OUT_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

STATS_DIR = os.path.join(OUT_DIR, "statistics")
os.makedirs(STATS_DIR, exist_ok=True)

LOG_FREQ = 50

print("=== RNAFOLD_DANEEL_PRO_ENHANCED (SMART PARSER + STATISTICS) ===")
print("[INFO] Input FASTA:", IN_FA)

if not os.path.exists(IN_FA):
    raise FileNotFoundError(f"Missing {IN_FA}")


# ================================================================
# USER INPUT (MFE RANGE)
# ================================================================
try:
    MIN_MFE = float(input("Enter MIN MFE (e.g., -60): ").strip())
    MAX_MFE = float(input("Enter MAX MFE (e.g., -15): ").strip())
except:
    raise ValueError("Invalid MFE values!")

print(f"[INFO] Filtering range: {MIN_MFE} .. {MAX_MFE}")


# ================================================================
# SMART HEADER PARSER
# ================================================================
def parse_header(header):
    """
    Smart FASTA header parser that tolerates ANY number of '|' or '='.
    Returns CoordID, length, gc_value, isoforms_list
    """

    parts = header.split("|")

    cid = parts[0]                      # always CoordID
    length = None
    gc_val = None
    iso_list = []

    for p in parts[1:]:
        if p.startswith("len="):
            length = p.split("=", 1)[1]
        elif p.startswith("gc="):
            gc_val = p.split("=", 1)[1]
        elif p.startswith("isoforms="):
            iso_list.append(p.split("=", 1)[1])
        else:
            # leftover fragments (caused by extra | in original header)
            iso_list.append(p)

    iso = ";".join(iso_list)
    return cid, length, gc_val, iso


# ================================================================
# ENHANCED STATISTICS FUNCTIONS
# ================================================================

def calculate_comprehensive_stats(results, in_range, out_range, MIN_MFE, MAX_MFE):
    """Kapsamlı MFE istatistiklerini hesapla"""
    stats = {}
    
    # Basic counts
    stats['total'] = len(results)
    stats['in_range'] = len(in_range)
    stats['out_range'] = len(out_range)
    stats['min_mfe'] = MIN_MFE
    stats['max_mfe'] = MAX_MFE
    
    # MFE values for all sequences
    mfe_vals = [r["MFE"] for r in results if r["MFE"] is not None]
    stats['all_mfe'] = mfe_vals
    
    if mfe_vals:
        stats['all_mean'] = statistics.mean(mfe_vals)
        stats['all_median'] = statistics.median(mfe_vals)
        stats['all_std'] = statistics.stdev(mfe_vals) if len(mfe_vals) > 1 else 0
        stats['all_min'] = min(mfe_vals)
        stats['all_max'] = max(mfe_vals)
        stats['all_percentiles'] = {
            p: np.percentile(mfe_vals, p) for p in [10, 25, 50, 75, 90, 95, 99]
        }
    else:
        stats.update({f'all_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['all_percentiles'] = {}
    
    # MFE values for in-range sequences
    in_range_mfe = [r["MFE"] for r in in_range if r["MFE"] is not None]
    stats['in_range_mfe'] = in_range_mfe
    
    if in_range_mfe:
        stats['in_range_mean'] = statistics.mean(in_range_mfe)
        stats['in_range_median'] = statistics.median(in_range_mfe)
        stats['in_range_std'] = statistics.stdev(in_range_mfe) if len(in_range_mfe) > 1 else 0
        stats['in_range_min'] = min(in_range_mfe)
        stats['in_range_max'] = max(in_range_mfe)
        stats['in_range_percentiles'] = {
            p: np.percentile(in_range_mfe, p) for p in [10, 25, 50, 75, 90, 95, 99]
        }
        
        # MFE distribution categories for in-range sequences
        stats['mfe_categories'] = {
            'very_stable': len([mfe for mfe in in_range_mfe if mfe < -40]),
            'stable': len([mfe for mfe in in_range_mfe if -40 <= mfe < -25]),
            'moderate': len([mfe for mfe in in_range_mfe if -25 <= mfe < -10]),
            'unstable': len([mfe for mfe in in_range_mfe if mfe >= -10])
        }
    else:
        stats.update({f'in_range_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['in_range_percentiles'] = {}
        stats['mfe_categories'] = {}
    
    # Filtering efficiency
    stats['retention_rate'] = (len(in_range) / len(results)) * 100 if len(results) > 0 else 0
    stats['filtering_rate'] = (len(out_range) / len(results)) * 100 if len(results) > 0 else 0
    
    # Structure statistics
    dot_lengths = [len(r["DotBracket"]) for r in results if r["DotBracket"]]
    if dot_lengths:
        stats['dot_mean_length'] = statistics.mean(dot_lengths)
        stats['dot_median_length'] = statistics.median(dot_lengths)
    
    # Sequence length vs MFE correlation
    lengths = [r["Length"] for r in results if r["MFE"] is not None]
    mfes = [r["MFE"] for r in results if r["MFE"] is not None]
    if lengths and mfes:
        stats['length_mfe_correlation'] = np.corrcoef(lengths, mfes)[0,1]
    else:
        stats['length_mfe_correlation'] = 0
    
    return stats

def generate_enhanced_plots(results, in_range, out_range, MIN_MFE, MAX_MFE):
    """Gelişmiş MFE grafikleri oluştur"""
    
    mfe_vals = [r["MFE"] for r in results if r["MFE"] is not None]
    in_range_mfe = [r["MFE"] for r in in_range if r["MFE"] is not None]
    
    # 1. Comparative MFE Histogram
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    if mfe_vals and in_range_mfe:
        plt.hist(mfe_vals, bins=40, alpha=0.7, color='lightblue', edgecolor='black', label='All Sequences')
        plt.hist(in_range_mfe, bins=40, alpha=0.7, color='lightgreen', edgecolor='black', label='In-Range Sequences')
        plt.axvline(MIN_MFE, color='red', linestyle='--', linewidth=2, label=f'Min ({MIN_MFE})')
        plt.axvline(MAX_MFE, color='red', linestyle='--', linewidth=2, label=f'Max ({MAX_MFE})')
        plt.xlabel('MFE (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('MFE Distribution: All vs In-Range Sequences')
        plt.legend()
    
    # 2. Category Distribution with percentages
    plt.subplot(2, 2, 2)
    categories = [f'Too Stable\n(<{MIN_MFE})', f'In Range\n({MIN_MFE}-{MAX_MFE})', f'Too Unstable\n(>{MAX_MFE})']
    counts = [len([r for r in out_range if r["MFE"] is not None and r["MFE"] < MIN_MFE]),
              len(in_range),
              len([r for r in out_range if r["MFE"] is not None and r["MFE"] > MAX_MFE])]
    percentages = [(c/len(results))*100 for c in counts]
    colors = ['blue', 'green', 'orange']
    
    bars = plt.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
    plt.ylabel('Count')
    plt.title('MFE Category Distribution')
    
    # Add value labels on bars
    for bar, count, percentage in zip(bars, counts, percentages):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01, 
                f'{count:,}\n({percentage:.1f}%)', 
                ha='center', va='bottom', fontsize=9)
    
    # 3. Boxplot comparison
    plt.subplot(2, 2, 3)
    if mfe_vals and in_range_mfe:
        data = [mfe_vals, in_range_mfe]
        labels = ['All Sequences', 'In-Range Sequences']
        plt.boxplot(data, labels=labels)
        plt.ylabel('MFE (kcal/mol)')
        plt.title('MFE Distribution Comparison')
    
    # 4. MFE distribution in in-range sequences
    plt.subplot(2, 2, 4)
    if in_range_mfe:
        plt.hist(in_range_mfe, bins=30, alpha=0.7, color='green', edgecolor='black')
        if in_range_mfe:
            plt.axvline(statistics.mean(in_range_mfe), color='red', linestyle='-', 
                       label=f'Mean: {statistics.mean(in_range_mfe):.1f}')
            plt.axvline(statistics.median(in_range_mfe), color='orange', linestyle='-', 
                       label=f'Median: {statistics.median(in_range_mfe):.1f}')
        plt.xlabel('MFE (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title('In-Range Sequences MFE Distribution')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "comprehensive_mfe_analysis.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Detailed in-range sequences MFE histogram
    plt.figure(figsize=(10, 6))
    if in_range_mfe:
        plt.hist(in_range_mfe, bins=40, alpha=0.7, color='darkgreen', edgecolor='black')
        plt.xlabel('MFE (kcal/mol)')
        plt.ylabel('Frequency')
        plt.title(f'Detailed MFE Distribution of In-Range Sequences ({MIN_MFE}-{MAX_MFE})')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "in_range_mfe_detail.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. Filtering efficiency pie chart
    plt.figure(figsize=(8, 6))
    labels = [f'In Range\n({len(in_range):,})', 
              f'Too Stable\n({len([r for r in out_range if r["MFE"] is not None and r["MFE"] < MIN_MFE]):,})', 
              f'Too Unstable\n({len([r for r in out_range if r["MFE"] is not None and r["MFE"] > MAX_MFE]):,})']
    sizes = [len(in_range),
             len([r for r in out_range if r["MFE"] is not None and r["MFE"] < MIN_MFE]),
             len([r for r in out_range if r["MFE"] is not None and r["MFE"] > MAX_MFE])]
    colors = ['green', 'blue', 'orange']
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    plt.title('MFE Filtering Efficiency')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "mfe_filtering_efficiency.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 7. Length vs MFE scatter plot
    plt.figure(figsize=(10, 6))
    lengths = [r["Length"] for r in results if r["MFE"] is not None]
    mfes = [r["MFE"] for r in results if r["MFE"] is not None]
    if lengths and mfes:
        plt.scatter(lengths, mfes, alpha=0.6, s=20)
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('MFE (kcal/mol)')
        plt.title(f'Sequence Length vs MFE (r = {np.corrcoef(lengths, mfes)[0,1]:.3f})')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "length_vs_mfe.png"), dpi=300, bbox_inches='tight')
        plt.close()

def save_detailed_statistics_report(stats, out_dir):
    """Detaylı MFE istatistik raporu oluştur"""
    report_path = os.path.join(out_dir, "detailed_mfe_statistics_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("RNAFOLD MFE FILTERING - DETAILED STATISTICAL REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        # Overview
        f.write("FILTERING OVERVIEW:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total Sequences Processed: {stats['total']:,}\n")
        f.write(f"Accepted Sequences ({stats['min_mfe']}-{stats['max_mfe']}): {stats['in_range']:,}\n")
        f.write(f"Filtered Out: {stats['out_range']:,}\n")
        f.write(f"  - Too Stable (<{stats['min_mfe']}): {len([mfe for mfe in stats['all_mfe'] if mfe < stats['min_mfe']]):,}\n")
        f.write(f"  - Too Unstable (>{stats['max_mfe']}): {len([mfe for mfe in stats['all_mfe'] if mfe > stats['max_mfe']]):,}\n")
        f.write(f"Retention Rate: {stats['retention_rate']:.2f}%\n")
        f.write(f"Filtering Rate: {stats['filtering_rate']:.2f}%\n\n")
        
        # All sequences MFE statistics
        f.write("ALL SEQUENCES MFE STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats.get('all_mean', 0) != 0:
            f.write(f"Mean MFE: {stats['all_mean']:.3f} kcal/mol\n")
            f.write(f"Median MFE: {stats['all_median']:.3f} kcal/mol\n")
            f.write(f"Standard Deviation: {stats['all_std']:.3f} kcal/mol\n")
            f.write(f"Range: {stats['all_min']:.3f} - {stats['all_max']:.3f} kcal/mol\n\n")
            
            f.write("PERCENTILES (All Sequences):\n")
            for p, value in stats['all_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.3f} kcal/mol\n")
            f.write("\n")
        else:
            f.write("No MFE data available\n\n")
        
        # In-range sequences MFE statistics
        f.write("IN-RANGE SEQUENCES MFE STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats.get('in_range_mean', 0) != 0:
            f.write(f"Mean MFE: {stats['in_range_mean']:.3f} kcal/mol\n")
            f.write(f"Median MFE: {stats['in_range_median']:.3f} kcal/mol\n")
            f.write(f"Standard Deviation: {stats['in_range_std']:.3f} kcal/mol\n")
            f.write(f"Range: {stats['in_range_min']:.3f} - {stats['in_range_max']:.3f} kcal/mol\n\n")
            
            f.write("PERCENTILES (In-Range Sequences):\n")
            for p, value in stats['in_range_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.3f} kcal/mol\n")
            f.write("\n")
            
            # MFE categories for in-range sequences
            if 'mfe_categories' in stats:
                f.write("MFE STABILITY CATEGORIES (In-Range Sequences):\n")
                for cat, count in stats['mfe_categories'].items():
                    percentage = (count / stats['in_range']) * 100 if stats['in_range'] > 0 else 0
                    f.write(f"  {cat.replace('_', ' ').title()}: {count:,} ({percentage:.1f}%)\n")
        else:
            f.write("No in-range sequences available\n\n")
        
        # Additional statistics
        f.write("ADDITIONAL STATISTICS:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Length-MFE Correlation: {stats['length_mfe_correlation']:.3f}\n")
        if 'dot_mean_length' in stats:
            f.write(f"Mean Dot-Bracket Length: {stats['dot_mean_length']:.1f}\n")
            f.write(f"Median Dot-Bracket Length: {stats['dot_median_length']:.1f}\n")
    
    print(f"[STATS] Detailed MFE statistics report saved to {report_path}")

# ================================================================
# LOAD FASTA
# ================================================================
def load_fasta(path):
    d = {}
    cur = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    d[cur] = "".join(buf)
                cur = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if cur:
            d[cur] = "".join(buf)
    return d

print("[INFO] Loading FASTA...")
fasta = load_fasta(IN_FA)
total = len(fasta)
print(f"[INFO] Loaded {total} records.")


# ================================================================
# RNAfold EXEC
# ================================================================
def run_rnafold(seq):
    p = subprocess.Popen(
        ["RNAfold", "--noPS"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    out, err = p.communicate(seq + "\n")
    lines = out.strip().split("\n")

    if len(lines) < 2:
        return "", None

    struct_line = lines[1]

    try:
        mfe_str = struct_line[struct_line.rfind("(")+1:struct_line.rfind(")")]
        mfe = float(mfe_str)
    except:
        mfe = None

    dot = struct_line.split(" ")[0]
    return dot, mfe


# ================================================================
# PROCESS SEQUENCES
# ================================================================
results = []
start = time.time()

print("\n[STEP] Running RNAfold (real-time)...\n")

for i, (header, seq) in enumerate(fasta.items(), 1):

    cid, L, gc_val, iso_str = parse_header(header)
    dot, mfe = run_rnafold(seq)

    if i % LOG_FREQ == 1:
        elapsed = time.time() - start
        pct = i / total * 100
        eps = i / elapsed if elapsed > 0 else 0
        eta = (total - i) / eps if eps > 0 else 0

        if eta < 60:
            eta_str = f"{eta:.0f}s"
        elif eta < 3600:
            eta_str = f"{eta/60:.1f}m"
        else:
            eta_str = f"{eta/3600:.1f}h"

        print(f"  → {i:,}/{total:,} ({pct:.1f}%) ETA {eta_str} - {cid}")

    results.append({
        "CoordID": cid,
        "Sequence": seq,
        "Length": len(seq),
        "DotBracket": dot,
        "MFE": mfe,
        "Isoforms": iso_str
    })

print("[INFO] RNAfold complete.")


# ================================================================
# FILTER BY MFE
# ================================================================
in_range  = []
out_range = []

for r in results:
    if r["MFE"] is None:
        out_range.append(r)
    elif MIN_MFE <= r["MFE"] <= MAX_MFE:
        in_range.append(r)
    else:
        out_range.append(r)


# ================================================================
# ENHANCED STATISTICS GENERATION
# ================================================================
print("\n" + "="*60)
print("GENERATING ENHANCED MFE STATISTICS")
print("="*60)

# Calculate enhanced MFE statistics
stats = calculate_comprehensive_stats(results, in_range, out_range, MIN_MFE, MAX_MFE)

# Generate enhanced MFE plots
generate_enhanced_plots(results, in_range, out_range, MIN_MFE, MAX_MFE)

# Save detailed MFE reports
save_detailed_statistics_report(stats, STATS_DIR)


# ================================================================
# WRITE TSV + FASTA 
# ================================================================
def write_tsv(path, rows):
    with open(path, "w", newline="") as f:
        fn = ["CoordID", "Sequence", "Length", "DotBracket", "MFE", "Isoforms"]
        w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

def write_fa(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write(
                f">{r['CoordID']}|len={r['Length']}|mfe={r['MFE']}|isoforms={r['Isoforms']}\n"
            )
            seq = r["Sequence"]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


# paths (ORIGINAL)
tsv_in  = os.path.join(OUT_DIR, "rna_inrange.tsv")
fa_in   = os.path.join(OUT_DIR, "rna_inrange.fa")
tsv_out = os.path.join(OUT_DIR, "rna_outrange.tsv")
fa_out  = os.path.join(OUT_DIR, "rna_outrange.fa")

write_tsv(tsv_in, in_range)
write_fa(fa_in, in_range)

write_tsv(tsv_out, out_range)
write_fa(fa_out, out_range)


# ================================================================
# ORIGINAL SUMMARY 
# ================================================================
mfe_vals = [r["MFE"] for r in results if r["MFE"] is not None]

mean_mfe = statistics.mean(mfe_vals)
med_mfe  = statistics.median(mfe_vals)
std_mfe  = statistics.pstdev(mfe_vals)

# Write summary file (ORIGINAL)
with open(os.path.join(OUT_DIR, "summary.tsv"), "w") as f:
    f.write("Metric\tValue\n")
    f.write(f"Total\t{len(results)}\n")
    f.write(f"InRange\t{len(in_range)}\n")
    f.write(f"OutRange\t{len(out_range)}\n")
    f.write(f"MeanMFE\t{mean_mfe:.3f}\n")
    f.write(f"MedianMFE\t{med_mfe:.3f}\n")
    f.write(f"StdMFE\t{std_mfe:.3f}\n")


# ================================================================
# ORIGINAL PLOTS
# ================================================================
plt.figure(figsize=(8,5))
plt.hist(mfe_vals, bins=40)
plt.title("MFE Histogram")
plt.xlabel("MFE")
plt.ylabel("Count")
plt.savefig(os.path.join(PLOT_DIR, "mfe_hist.png"), dpi=200)
plt.close()

plt.figure(figsize=(8,5))
plt.boxplot(mfe_vals)
plt.title("MFE Boxplot")
plt.ylabel("MFE")
plt.savefig(os.path.join(PLOT_DIR, "mfe_box.png"), dpi=200)
plt.close()


# ================================================================
# ENHANCED TERMINAL SUMMARY
# ================================================================
print("\n" + "="*65)
print("RNAFOLD SUMMARY - ENHANCED")
print("="*65)
print(f"Total sequences       : {len(results):,}")
print(f"In-range MFE          : {len(in_range):,} ({stats['retention_rate']:.1f}%)")
print(f"Out-of-range MFE      : {len(out_range):,} ({stats['filtering_rate']:.1f}%)")
print(f"  - Too stable        : {len([r for r in out_range if r['MFE'] is not None and r['MFE'] < MIN_MFE]):,}")
print(f"  - Too unstable      : {len([r for r in out_range if r['MFE'] is not None and r['MFE'] > MAX_MFE]):,}")
print()
print("MFE STATISTICS (In-Range Sequences):")
print(f"  Mean MFE            : {stats['in_range_mean']:.3f} kcal/mol")
print(f"  Median MFE          : {stats['in_range_median']:.3f} kcal/mol")
print(f"  Std deviation       : {stats['in_range_std']:.3f} kcal/mol")
print(f"  Range               : {stats['in_range_min']:.3f} - {stats['in_range_max']:.3f} kcal/mol")
print()
print("MFE STATISTICS (All Sequences):")
print(f"  Mean MFE            : {stats['all_mean']:.3f} kcal/mol")
print(f"  Median MFE          : {stats['all_median']:.3f} kcal/mol")
print(f"  Length-MFE Corr     : {stats['length_mfe_correlation']:.3f}")
print("="*65)

print(f"\n[✓] DONE — All RNAfold outputs written to {OUT_DIR}/")
print(f"[📊] Enhanced MFE statistics: {STATS_DIR}/")
print(f"[📈] Enhanced MFE plots: {PLOT_DIR}/")
print(f"[📋] Original outputs preserved")