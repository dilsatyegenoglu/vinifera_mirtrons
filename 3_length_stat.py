#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LENGTH_ENHANCED (Pipeline-Safe + Statistics)
Fully corrected version with enhanced statistics:
- Real-time filtering with ETA
- IN/OUT FASTA + TSV (pipeline-safe headers)
- Comprehensive statistical analysis
- Enhanced visualizations
- Detailed reports
"""

import os
import csv
import time
import statistics
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# ============================================================
# PATHS
# ============================================================

BASE = os.path.abspath(os.path.dirname(__file__))
IN_DIR = os.path.join(BASE, "consolidated_introns")
OUT_DIR = os.path.join(BASE, "length_filter_results")
PLOT_DIR = os.path.join(OUT_DIR, "plots")
STATS_DIR = os.path.join(OUT_DIR, "statistics")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(STATS_DIR, exist_ok=True)

TSV_IN = os.path.join(IN_DIR, "consolidated.tsv")

# ============================================================
# LOAD INPUT
# ============================================================

def load_rows():
    rows = []
    with open(TSV_IN) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

# ============================================================
# ENHANCED STATISTICS FUNCTIONS
# ============================================================

def calculate_comprehensive_stats(all_lengths, kept_lengths, too_short, too_long, kept, total, MIN_L, MAX_L):
    """Calculation of Comprehensive Statistics"""
    stats = {}
    
    # Basic counts
    stats['total'] = total
    stats['kept'] = kept
    stats['filtered_out'] = total - kept
    stats['too_short'] = too_short
    stats['too_long'] = too_long
    stats['min_accepted'] = MIN_L
    stats['max_accepted'] = MAX_L
    
    # Length statistics for all introns
    if all_lengths:
        stats['all_mean'] = statistics.mean(all_lengths)
        stats['all_median'] = statistics.median(all_lengths)
        stats['all_std'] = statistics.stdev(all_lengths) if len(all_lengths) > 1 else 0
        stats['all_min'] = min(all_lengths)
        stats['all_max'] = max(all_lengths)
        stats['all_percentiles'] = {
            p: np.percentile(all_lengths, p) for p in [25, 50, 75, 90, 95, 99]
        }
    else:
        stats.update({f'all_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['all_percentiles'] = {}
    
    # Length statistics for kept introns
    if kept_lengths:
        stats['kept_mean'] = statistics.mean(kept_lengths)
        stats['kept_median'] = statistics.median(kept_lengths)
        stats['kept_std'] = statistics.stdev(kept_lengths) if len(kept_lengths) > 1 else 0
        stats['kept_min'] = min(kept_lengths)
        stats['kept_max'] = max(kept_lengths)
        stats['kept_percentiles'] = {
            p: np.percentile(kept_lengths, p) for p in [25, 50, 75, 90, 95, 99]
        }
    else:
        stats.update({f'kept_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['kept_percentiles'] = {}
    
    # Length distribution categories
    if kept_lengths:
        stats['length_categories'] = {
            'very_short': len([l for l in kept_lengths if l < 50]),
            'short': len([l for l in kept_lengths if 50 <= l < 100]),
            'medium': len([l for l in kept_lengths if 100 <= l < 500]),
            'long': len([l for l in kept_lengths if 500 <= l < 2000]),
            'very_long': len([l for l in kept_lengths if l >= 2000])
        }
    
    # Filtering efficiency
    stats['retention_rate'] = (kept / total) * 100 if total > 0 else 0
    stats['filtering_rate'] = ((total - kept) / total) * 100 if total > 0 else 0
    
    return stats

def generate_enhanced_plots(all_lengths, kept_lengths, too_short, too_long, kept, total, MIN_L, MAX_L):
    """Gelişmiş grafikler oluştur"""
    
    # 1. Comparative Histogram
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    if all_lengths:
        plt.hist(all_lengths, bins=50, alpha=0.7, color='lightblue', edgecolor='black', label='All Introns')
        plt.hist(kept_lengths, bins=50, alpha=0.7, color='lightgreen', edgecolor='black', label='Kept Introns')
        plt.axvline(MIN_L, color='red', linestyle='--', linewidth=2, label=f'Min ({MIN_L})')
        plt.axvline(MAX_L, color='red', linestyle='--', linewidth=2, label=f'Max ({MAX_L})')
        plt.xlabel('Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Length Distribution: All vs Kept Introns')
        plt.legend()
        plt.yscale('log')
    
    # 2. Category Distribution with percentages
    plt.subplot(2, 2, 2)
    categories = [f'Too Short\n(<{MIN_L})', f'In Range\n({MIN_L}-{MAX_L})', f'Too Long\n(>{MAX_L})']
    counts = [too_short, kept, too_long]
    percentages = [(c/total)*100 for c in counts]
    colors = ['red', 'green', 'orange']
    
    bars = plt.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
    plt.ylabel('Count')
    plt.title('Length Category Distribution')
    
    # Add value labels on bars
    for bar, count, percentage in zip(bars, counts, percentages):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01, 
                f'{count:,}\n({percentage:.1f}%)', 
                ha='center', va='bottom', fontsize=9)
    
    # 3. Boxplot comparison
    plt.subplot(2, 2, 3)
    if all_lengths and kept_lengths:
        data = [all_lengths, kept_lengths]
        labels = ['All Introns', 'Kept Introns']
        plt.boxplot(data, labels=labels)
        plt.ylabel('Length (bp)')
        plt.title('Length Distribution Comparison')
    
    # 4. Length distribution in kept introns
    plt.subplot(2, 2, 4)
    if kept_lengths:
        plt.hist(kept_lengths, bins=30, alpha=0.7, color='green', edgecolor='black')
        if kept_lengths:
            plt.axvline(statistics.mean(kept_lengths), color='red', linestyle='-', 
                       label=f'Mean: {statistics.mean(kept_lengths):.1f}')
            plt.axvline(statistics.median(kept_lengths), color='orange', linestyle='-', 
                       label=f'Median: {statistics.median(kept_lengths):.1f}')
        plt.xlabel('Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Kept Introns Length Distribution')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "comprehensive_analysis.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Detailed kept introns histogram
    plt.figure(figsize=(10, 6))
    if kept_lengths:
        plt.hist(kept_lengths, bins=50, alpha=0.7, color='darkgreen', edgecolor='black')
        plt.xlabel('Length (bp)')
        plt.ylabel('Frequency')
        plt.title(f'Detailed Length Distribution of Kept Introns ({MIN_L}-{MAX_L} bp)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "kept_introns_detail.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. Filtering efficiency pie chart
    plt.figure(figsize=(8, 6))
    labels = [f'Kept\n({kept:,})', f'Too Short\n({too_short:,})', f'Too Long\n({too_long:,})']
    sizes = [kept, too_short, too_long]
    colors = ['green', 'red', 'orange']
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    plt.title('Filtering Efficiency')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "filtering_efficiency.png"), dpi=300, bbox_inches='tight')
    plt.close()

def save_detailed_statistics_report(stats, out_dir):
    """Detailed statistics report"""
    report_path = os.path.join(out_dir, "detailed_statistics_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("LENGTH FILTERING - DETAILED STATISTICAL REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        # Overview
        f.write("FILTERING OVERVIEW:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total Introns Processed: {stats['total']:,}\n")
        f.write(f"Accepted Introns ({stats['min_accepted']}-{stats['max_accepted']} bp): {stats['kept']:,}\n")
        f.write(f"Filtered Out: {stats['filtered_out']:,}\n")
        f.write(f"  - Too Short (<{stats['min_accepted']} bp): {stats['too_short']:,}\n")
        f.write(f"  - Too Long (>{stats['max_accepted']} bp): {stats['too_long']:,}\n")
        f.write(f"Retention Rate: {stats['retention_rate']:.2f}%\n")
        f.write(f"Filtering Rate: {stats['filtering_rate']:.2f}%\n\n")
        
        # All introns statistics
        f.write("ALL INTRONS LENGTH STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats['all_lengths']:
            f.write(f"Mean Length: {stats['all_mean']:.2f} bp\n")
            f.write(f"Median Length: {stats['all_median']:.2f} bp\n")
            f.write(f"Standard Deviation: {stats['all_std']:.2f} bp\n")
            f.write(f"Range: {stats['all_min']:,} - {stats['all_max']:,} bp\n\n")
            
            f.write("PERCENTILES (All Introns):\n")
            for p, value in stats['all_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.2f} bp\n")
            f.write("\n")
        else:
            f.write("No length data available\n\n")
        
        # Kept introns statistics
        f.write("KEPT INTRONS LENGTH STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats['kept_lengths']:
            f.write(f"Mean Length: {stats['kept_mean']:.2f} bp\n")
            f.write(f"Median Length: {stats['kept_median']:.2f} bp\n")
            f.write(f"Standard Deviation: {stats['kept_std']:.2f} bp\n")
            f.write(f"Range: {stats['kept_min']:,} - {stats['kept_max']:,} bp\n\n")
            
            f.write("PERCENTILES (Kept Introns):\n")
            for p, value in stats['kept_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.2f} bp\n")
            f.write("\n")
            
            # Length categories for kept introns
            if 'length_categories' in stats:
                f.write("LENGTH CATEGORIES (Kept Introns):\n")
                for cat, count in stats['length_categories'].items():
                    percentage = (count / stats['kept']) * 100
                    f.write(f"  {cat.replace('_', ' ').title()}: {count:,} ({percentage:.1f}%)\n")
        else:
            f.write("No kept introns available\n")
    
    print(f"[STATS] Detailed statistics report saved to {report_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    print("=== LENGTH6_DANEEL_FINAL_ENHANCED (Real-time + Enhanced Stats) ===")

    rows = load_rows()
    total = len(rows)
    print(f"[INFO] Loaded {total} introns")

    MIN_L = int(input("Enter MIN length: ").strip())
    MAX_L = int(input("Enter MAX length: ").strip())

    # OUTPUT FILES (ORIJINAL YAPISI KORUNDU)
    in_tsv  = open(os.path.join(OUT_DIR, "inrange.tsv"), "w")
    in_fa   = open(os.path.join(OUT_DIR, "inrange.fa"), "w")
    out_tsv = open(os.path.join(OUT_DIR, "outrange.tsv"), "w")
    out_fa  = open(os.path.join(OUT_DIR, "outrange.fa"), "w")

    header = ["CoordID","Length","IsoformCount","IsoformList","Sequence","Status"]
    in_tsv.write("\t".join(header) + "\n")
    out_tsv.write("\t".join(header) + "\n")

    kept_lengths = []
    all_lengths = []
    too_short = 0
    too_long = 0
    kept = 0

    start = time.time()

    # ======================
    # REAL-TIME LOOP (ORIJINAL KOD)
    # ======================
    for i, r in enumerate(rows, 1):

        cid = r["CoordID"].strip()
        seq = r["Sequence"].strip()
        iso_list = r["IsoformList"]
        iso_count = r["IsoformCount"]

        try:
            L = len(seq)
        except:
            L = 0

        all_lengths.append(L)

        # ETA log
        if i % 500 == 1:
            elapsed = time.time() - start
            pct = (i/total)*100
            ips = i/elapsed if elapsed > 0 else 0
            rem = (total - i) / ips if ips > 0 else 0

            if rem < 60:
                eta = f"{rem:.0f}s"
            elif rem < 3600:
                eta = f"{rem/60:.1f}m"
            else:
                eta = f"{rem/3600:.1f}h"

            print(f"→ {i:,}/{total:,} ({pct:.1f}%)  ETA: {eta}  {cid}")

        # Filtering (ORIJINAL KOD)
        if L < MIN_L:
            too_short += 1
            status = "TOO_SHORT"
            out_tsv.write(f"{cid}\t{L}\t{iso_count}\t{iso_list}\t{seq}\t{status}\n")
            out_fa.write(f">{cid}|len={L}|isoforms={iso_list}\n")
            for j in range(0, len(seq), 80):
                out_fa.write(seq[j:j+80] + "\n")

        elif L > MAX_L:
            too_long += 1
            status = "TOO_LONG"
            out_tsv.write(f"{cid}\t{L}\t{iso_count}\t{iso_list}\t{seq}\t{status}\n")
            out_fa.write(f">{cid}|len={L}|isoforms={iso_list}\n")
            for j in range(0, len(seq), 80):
                out_fa.write(seq[j:j+80] + "\n")

        else:
            kept += 1
            kept_lengths.append(L)
            status = "IN_RANGE"
            in_tsv.write(f"{cid}\t{L}\t{iso_count}\t{iso_list}\t{seq}\t{status}\n")
            in_fa.write(f">{cid}|len={L}|isoforms={iso_list}\n")
            for j in range(0, len(seq), 80):
                in_fa.write(seq[j:j+80] + "\n")

    in_tsv.close()
    in_fa.close()
    out_tsv.close()
    out_fa.close()

    # ======================================================
    # ENHANCED STATISTICS
    # ======================================================
    print("\n" + "="*60)
    print("GENERATING ENHANCED STATISTICS")
    print("="*60)
    
    # Calculate comprehensive statistics
    stats = calculate_comprehensive_stats(all_lengths, kept_lengths, too_short, too_long, kept, total, MIN_L, MAX_L)
    stats['all_lengths'] = all_lengths
    stats['kept_lengths'] = kept_lengths
    
    # Generate enhanced plots
    generate_enhanced_plots(all_lengths, kept_lengths, too_short, too_long, kept, total, MIN_L, MAX_L)
    
    # Make detailed report
    save_detailed_statistics_report(stats, STATS_DIR)

    # ======================================================
    # ORIGINAL STATISTICS
    # ======================================================
    meanL = statistics.mean(kept_lengths) if kept_lengths else 0
    medL  = statistics.median(kept_lengths) if kept_lengths else 0
    stdL  = statistics.stdev(kept_lengths) if len(kept_lengths)>1 else 0

    # summary.tsv (ORIJINAL)
    with open(os.path.join(OUT_DIR, "summary.tsv"), "w") as f:
        f.write("Metric\tValue\n")
        f.write(f"Total\t{total}\n")
        f.write(f"MinAccepted\t{MIN_L}\n")
        f.write(f"MaxAccepted\t{MAX_L}\n")
        f.write(f"Kept\t{kept}\n")
        f.write(f"FilteredOut\t{total-kept}\n")
        f.write(f"TooShort\t{too_short}\n")
        f.write(f"TooLong\t{too_long}\n")
        f.write(f"MeanLength\t{meanL:.3f}\n")
        f.write(f"MedianLength\t{medL}\n")
        f.write(f"StdDevLength\t{stdL:.3f}\n")

    # pretty summary (ORIGINAL)
    with open(os.path.join(OUT_DIR, "summary_pretty.txt"), "w") as f:
        f.write("=============== LENGTH FILTER SUMMARY ===============\n")
        f.write(f"Total introns        : {total}\n")
        f.write(f"Accepted ({MIN_L}-{MAX_L}) : {kept}\n")
        f.write(f"Filtered out         : {total-kept}\n")
        f.write(f"   - Too short       : {too_short}\n")
        f.write(f"   - Too long        : {too_long}\n")
        f.write("\n")
        f.write(f"Mean length          : {meanL:.3f}\n")
        f.write(f"Median length        : {medL}\n")
        f.write(f"Std deviation        : {stdL:.3f}\n")
        f.write("=====================================================\n")

    # ======================================================
    # ORIGINAL PLOTS
    # ======================================================

    # Histogram
    plt.figure(figsize=(8,5))
    plt.hist(all_lengths, bins=50)
    plt.title("Length Distribution (All Introns)")
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "histogram.png"))
    plt.close()

    # Boxplot
    plt.figure(figsize=(5,5))
    plt.boxplot(all_lengths)
    plt.title("Length Boxplot")
    plt.ylabel("Length")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "boxplot.png"))
    plt.close()

    # Barplot
    plt.figure(figsize=(6,5))
    plt.bar(["Too Short","In Range","Too Long"], [too_short, kept, too_long])
    plt.title("Length Category Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "barplot.png"))
    plt.close()

    # Scatter
    plt.figure(figsize=(8,5))
    plt.scatter(range(len(all_lengths)), all_lengths, s=5)
    plt.title("Index vs Length")
    plt.xlabel("Index")
    plt.ylabel("Length")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "scatter.png"))
    plt.close()

    # ======================================================
    # ENHANCED TERMINAL SUMMARY
    # ======================================================

    print("\n" + "="*65)
    print("LENGTH FILTER SUMMARY - ENHANCED")
    print("="*65)
    print(f"Total introns          : {total:,}")
    print(f"Accepted ({MIN_L}-{MAX_L})     : {kept:,} ({stats['retention_rate']:.1f}%)")
    print(f"Filtered out           : {total-kept:,} ({stats['filtering_rate']:.1f}%)")
    print(f"   - Too short         : {too_short:,} ({(too_short/total)*100:.1f}%)")
    print(f"   - Too long          : {too_long:,} ({(too_long/total)*100:.1f}%)")
    print()
    print("LENGTH STATISTICS (Kept Introns):")
    print(f"  Mean length          : {stats['kept_mean']:.1f} bp")
    print(f"  Median length        : {stats['kept_median']:.1f} bp")
    print(f"  Std deviation        : {stats['kept_std']:.1f} bp")
    print(f"  Range                : {stats['kept_min']:,} - {stats['kept_max']:,} bp")
    print()
    print("LENGTH STATISTICS (All Introns):")
    print(f"  Mean length          : {stats['all_mean']:.1f} bp")
    print(f"  Median length        : {stats['all_median']:.1f} bp")
    print("="*65)

    print(f"\n[✓] DONE — All outputs written to {OUT_DIR}/")
    print(f"[📊] Enhanced statistics: {STATS_DIR}/")
    print(f"[📈] Enhanced plots: {PLOT_DIR}/")
    print(f"[📋] Original outputs preserved")

# ============================================================

if __name__ == "__main__":
    main()