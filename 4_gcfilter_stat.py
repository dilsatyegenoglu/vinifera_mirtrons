#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GC_FILTER_ENHANCED
Enhanced pipeline-safe GC filtering module with comprehensive statistics.
Input  : length_filter_results/inrange.tsv + inrange.fa
Output : gc_filter_results/  (tsv + fasta + plots + summaries + statistics)
Features:
 - Real-time progress with ETA
 - Min/Max GC threshold (user input)
 - IN/OUT TSV + FASTA (pipeline-safe FASTA headers)
 - Comprehensive statistical analysis
 - Enhanced visualizations
 - Detailed reports
 - All original outputs preserved
"""

import os
import csv
import time
import math
import statistics
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# ============================================================
# PATHS
# ============================================================

BASE = os.path.abspath(os.path.dirname(__file__))
IN_DIR  = os.path.join(BASE, "length_filter_results")
OUT_DIR = os.path.join(BASE, "gc_filter_results")
PLOT_DIR = os.path.join(OUT_DIR, "plots")
STATS_DIR = os.path.join(OUT_DIR, "statistics")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(STATS_DIR, exist_ok=True)

TSV_IN = os.path.join(IN_DIR, "inrange.tsv")
FA_IN  = os.path.join(IN_DIR, "inrange.fa")


# ============================================================
# LOAD TSV
# ============================================================

def load_rows():
    rows = []
    with open(TSV_IN) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows


# ============================================================
# FASTA READER
# ============================================================

def load_fasta(path):
    d = {}
    cur_id = None
    seq_buf = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    d[cur_id] = "".join(seq_buf)
                cur_id = line[1:].strip()  # full header
                seq_buf = []
            else:
                seq_buf.append(line.strip())

        if cur_id:
            d[cur_id] = "".join(seq_buf)

    return d


# ============================================================
# GC CALCULATOR
# ============================================================

def gc_percent(seq):
    if not seq:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return ((g + c) / len(seq)) * 100.0

# ============================================================
# ENHANCED STATISTICS FUNCTIONS
# ============================================================

def calculate_comprehensive_stats(all_gc, kept_gc, low_gc, high_gc, kept, total, MIN_GC, MAX_GC):
    """Kapsamlı GC istatistiklerini hesapla"""
    stats = {}
    
    # Basic counts
    stats['total'] = total
    stats['kept'] = kept
    stats['filtered_out'] = total - kept
    stats['low_gc'] = low_gc
    stats['high_gc'] = high_gc
    stats['min_gc'] = MIN_GC
    stats['max_gc'] = MAX_GC
    
    # GC statistics for all introns
    if all_gc:
        stats['all_mean'] = statistics.mean(all_gc)
        stats['all_median'] = statistics.median(all_gc)
        stats['all_std'] = statistics.stdev(all_gc) if len(all_gc) > 1 else 0
        stats['all_min'] = min(all_gc)
        stats['all_max'] = max(all_gc)
        stats['all_percentiles'] = {
            p: np.percentile(all_gc, p) for p in [10, 25, 50, 75, 90, 95, 99]
        }
    else:
        stats.update({f'all_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['all_percentiles'] = {}
    
    # GC statistics for kept introns
    if kept_gc:
        stats['kept_mean'] = statistics.mean(kept_gc)
        stats['kept_median'] = statistics.median(kept_gc)
        stats['kept_std'] = statistics.stdev(kept_gc) if len(kept_gc) > 1 else 0
        stats['kept_min'] = min(kept_gc)
        stats['kept_max'] = max(kept_gc)
        stats['kept_percentiles'] = {
            p: np.percentile(kept_gc, p) for p in [10, 25, 50, 75, 90, 95, 99]
        }
        
        # GC distribution categories for kept introns
        stats['gc_categories'] = {
            'very_low': len([gc for gc in kept_gc if gc < 20]),
            'low': len([gc for gc in kept_gc if 20 <= gc < 35]),
            'medium': len([gc for gc in kept_gc if 35 <= gc < 50]),
            'high': len([gc for gc in kept_gc if 50 <= gc < 65]),
            'very_high': len([gc for gc in kept_gc if gc >= 65])
        }
    else:
        stats.update({f'kept_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['kept_percentiles'] = {}
        stats['gc_categories'] = {}
    
    # Filtering efficiency
    stats['retention_rate'] = (kept / total) * 100 if total > 0 else 0
    stats['filtering_rate'] = ((total - kept) / total) * 100 if total > 0 else 0
    stats['low_gc_rate'] = (low_gc / total) * 100 if total > 0 else 0
    stats['high_gc_rate'] = (high_gc / total) * 100 if total > 0 else 0
    
    return stats

def generate_enhanced_plots(all_gc, kept_gc, low_gc, high_gc, kept, total, MIN_GC, MAX_GC):
    """Gelişmiş GC grafikleri oluştur"""
    
    # 1. Comparative GC Histogram
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    if all_gc:
        plt.hist(all_gc, bins=40, alpha=0.7, color='lightblue', edgecolor='black', label='All Introns')
        plt.hist(kept_gc, bins=40, alpha=0.7, color='lightgreen', edgecolor='black', label='Kept Introns')
        plt.axvline(MIN_GC, color='red', linestyle='--', linewidth=2, label=f'Min ({MIN_GC}%)')
        plt.axvline(MAX_GC, color='red', linestyle='--', linewidth=2, label=f'Max ({MAX_GC}%)')
        plt.xlabel('GC Content (%)')
        plt.ylabel('Frequency')
        plt.title('GC Distribution: All vs Kept Introns')
        plt.legend()
    
    # 2. Category Distribution with percentages
    plt.subplot(2, 2, 2)
    categories = [f'Too Low\n(<{MIN_GC}%)', f'In Range\n({MIN_GC}-{MAX_GC}%)', f'Too High\n(>{MAX_GC}%)']
    counts = [low_gc, kept, high_gc]
    percentages = [(c/total)*100 for c in counts]
    colors = ['red', 'green', 'orange']
    
    bars = plt.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
    plt.ylabel('Count')
    plt.title('GC Category Distribution')
    
    # Add value labels on bars
    for bar, count, percentage in zip(bars, counts, percentages):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01, 
                f'{count:,}\n({percentage:.1f}%)', 
                ha='center', va='bottom', fontsize=9)
    
    # 3. Boxplot comparison
    plt.subplot(2, 2, 3)
    if all_gc and kept_gc:
        data = [all_gc, kept_gc]
        labels = ['All Introns', 'Kept Introns']
        plt.boxplot(data, labels=labels)
        plt.ylabel('GC Content (%)')
        plt.title('GC Distribution Comparison')
    
    # 4. GC distribution in kept introns
    plt.subplot(2, 2, 4)
    if kept_gc:
        plt.hist(kept_gc, bins=30, alpha=0.7, color='green', edgecolor='black')
        if kept_gc:
            plt.axvline(statistics.mean(kept_gc), color='red', linestyle='-', 
                       label=f'Mean: {statistics.mean(kept_gc):.1f}%')
            plt.axvline(statistics.median(kept_gc), color='orange', linestyle='-', 
                       label=f'Median: {statistics.median(kept_gc):.1f}%')
        plt.xlabel('GC Content (%)')
        plt.ylabel('Frequency')
        plt.title('Kept Introns GC Distribution')
        plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "comprehensive_gc_analysis.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 5. Detailed kept introns GC histogram
    plt.figure(figsize=(10, 6))
    if kept_gc:
        plt.hist(kept_gc, bins=40, alpha=0.7, color='darkgreen', edgecolor='black')
        plt.xlabel('GC Content (%)')
        plt.ylabel('Frequency')
        plt.title(f'Detailed GC Distribution of Kept Introns ({MIN_GC}-{MAX_GC}%)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "kept_introns_gc_detail.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # 6. Filtering efficiency pie chart
    plt.figure(figsize=(8, 6))
    labels = [f'Kept\n({kept:,})', f'Too Low GC\n({low_gc:,})', f'Too High GC\n({high_gc:,})']
    sizes = [kept, low_gc, high_gc]
    colors = ['green', 'red', 'orange']
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    plt.title('GC Filtering Efficiency')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "gc_filtering_efficiency.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 7. GC distribution by category (kept introns only)
    plt.figure(figsize=(10, 6))
    if kept_gc:
        gc_categories = {
            'Very Low (<20%)': len([gc for gc in kept_gc if gc < 20]),
            'Low (20-35%)': len([gc for gc in kept_gc if 20 <= gc < 35]),
            'Medium (35-50%)': len([gc for gc in kept_gc if 35 <= gc < 50]),
            'High (50-65%)': len([gc for gc in kept_gc if 50 <= gc < 65]),
            'Very High (≥65%)': len([gc for gc in kept_gc if gc >= 65])
        }
        
        categories, counts = zip(*gc_categories.items())
        plt.bar(categories, counts, color='purple', alpha=0.7, edgecolor='black')
        plt.xlabel('GC Content Category')
        plt.ylabel('Count')
        plt.title('GC Content Distribution in Kept Introns')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(PLOT_DIR, "gc_categories_kept.png"), dpi=300, bbox_inches='tight')
        plt.close()

def save_detailed_statistics_report(stats, out_dir):
    """Detaylı GC istatistik raporu oluştur"""
    report_path = os.path.join(out_dir, "detailed_gc_statistics_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("GC FILTERING - DETAILED STATISTICAL REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        # Overview
        f.write("FILTERING OVERVIEW:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total Introns Processed: {stats['total']:,}\n")
        f.write(f"Accepted Introns ({stats['min_gc']}-{stats['max_gc']}%): {stats['kept']:,}\n")
        f.write(f"Filtered Out: {stats['filtered_out']:,}\n")
        f.write(f"  - Too Low GC (<{stats['min_gc']}%): {stats['low_gc']:,}\n")
        f.write(f"  - Too High GC (>{stats['max_gc']}%): {stats['high_gc']:,}\n")
        f.write(f"Retention Rate: {stats['retention_rate']:.2f}%\n")
        f.write(f"Filtering Rate: {stats['filtering_rate']:.2f}%\n")
        f.write(f"Low GC Rate: {stats['low_gc_rate']:.2f}%\n")
        f.write(f"High GC Rate: {stats['high_gc_rate']:.2f}%\n\n")
        
        # All introns GC statistics
        f.write("ALL INTRONS GC STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats.get('all_mean', 0) > 0:
            f.write(f"Mean GC Content: {stats['all_mean']:.2f}%\n")
            f.write(f"Median GC Content: {stats['all_median']:.2f}%\n")
            f.write(f"Standard Deviation: {stats['all_std']:.2f}%\n")
            f.write(f"Range: {stats['all_min']:.2f}% - {stats['all_max']:.2f}%\n\n")
            
            f.write("PERCENTILES (All Introns):\n")
            for p, value in stats['all_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.2f}%\n")
            f.write("\n")
        else:
            f.write("No GC data available\n\n")
        
        # Kept introns GC statistics
        f.write("KEPT INTRONS GC STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats.get('kept_mean', 0) > 0:
            f.write(f"Mean GC Content: {stats['kept_mean']:.2f}%\n")
            f.write(f"Median GC Content: {stats['kept_median']:.2f}%\n")
            f.write(f"Standard Deviation: {stats['kept_std']:.2f}%\n")
            f.write(f"Range: {stats['kept_min']:.2f}% - {stats['kept_max']:.2f}%\n\n")
            
            f.write("PERCENTILES (Kept Introns):\n")
            for p, value in stats['kept_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.2f}%\n")
            f.write("\n")
            
            # GC categories for kept introns
            if 'gc_categories' in stats:
                f.write("GC CONTENT CATEGORIES (Kept Introns):\n")
                for cat, count in stats['gc_categories'].items():
                    percentage = (count / stats['kept']) * 100 if stats['kept'] > 0 else 0
                    f.write(f"  {cat.replace('_', ' ').title()}: {count:,} ({percentage:.1f}%)\n")
        else:
            f.write("No kept introns available\n")
    
    print(f"[STATS] Detailed GC statistics report saved to {report_path}")

# ============================================================
# MAIN
# ============================================================

def main():
    print("=== GC_FILTER_ENHANCED (Real-time + Enhanced Stats) ===")

    rows = load_rows()
    total = len(rows)
    print(f"[INFO] Loaded {total} introns from inrange.tsv")

    fasta_dict = load_fasta(FA_IN)

    MIN_GC = float(input("Enter MIN GC%: ").strip())
    MAX_GC = float(input("Enter MAX GC%: ").strip())

    # OUTPUT FILES 
    in_tsv  = open(os.path.join(OUT_DIR, "gc_inrange.tsv"), "w")
    in_fa   = open(os.path.join(OUT_DIR, "gc_inrange.fa"), "w")
    out_tsv = open(os.path.join(OUT_DIR, "gc_outrange.tsv"), "w")
    out_fa  = open(os.path.join(OUT_DIR, "gc_outrange.fa"), "w")

    header = ["CoordID","Length","GC_Percent","IsoformCount","IsoformList","Sequence","Status"]
    in_tsv.write("\t".join(header) + "\n")
    out_tsv.write("\t".join(header) + "\n")

    kept_gc = []
    all_gc = []
    low_gc = 0
    high_gc = 0
    kept = 0

    start = time.time()

    # REAL-TIME LOOP
    for i, r in enumerate(rows, 1):

        cid = r["CoordID"]
        seq = r["Sequence"]
        iso_list = r["IsoformList"]
        iso_count = r["IsoformCount"]
        length = len(seq)

        # GC%
        gcp = gc_percent(seq)
        all_gc.append(gcp)

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

        # FILTERING
        if gcp < MIN_GC:
            low_gc += 1
            status = "GC_TOO_LOW"
            out_tsv.write(f"{cid}\t{length}\t{gcp:.3f}\t{iso_count}\t{iso_list}\t{seq}\t{status}\n")
            out_fa.write(f">{cid}|len={length}|gc={gcp:.3f}|isoforms={iso_list}\n")
            for j in range(0, len(seq), 80):
                out_fa.write(seq[j:j+80] + "\n")

        elif gcp > MAX_GC:
            high_gc += 1
            status = "GC_TOO_HIGH"
            out_tsv.write(f"{cid}\t{length}\t{gcp:.3f}\t{iso_count}\t{iso_list}\t{seq}\t{status}\n")
            out_fa.write(f">{cid}|len={length}|gc={gcp:.3f}|isoforms={iso_list}\n")
            for j in range(0, len(seq), 80):
                out_fa.write(seq[j:j+80] + "\n")

        else:
            kept += 1
            kept_gc.append(gcp)
            status = "IN_RANGE"
            in_tsv.write(f"{cid}\t{length}\t{gcp:.3f}\t{iso_count}\t{iso_list}\t{seq}\t{status}\n")
            in_fa.write(f">{cid}|len={length}|gc={gcp:.3f}|isoforms={iso_list}\n")
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
    print("GENERATING ENHANCED GC STATISTICS")
    print("="*60)
    
    # Calculate comprehensive GC statistics
    stats = calculate_comprehensive_stats(all_gc, kept_gc, low_gc, high_gc, kept, total, MIN_GC, MAX_GC)
    
    # Generate enhanced GC plots
    generate_enhanced_plots(all_gc, kept_gc, low_gc, high_gc, kept, total, MIN_GC, MAX_GC)
    
    # Save detailed GC statistics
    save_detailed_statistics_report(stats, STATS_DIR)

    # ======================================================
    # ORIGINAL STATISTICS
    # ======================================================
    mean_gc = statistics.mean(kept_gc) if kept_gc else 0
    med_gc  = statistics.median(kept_gc) if kept_gc else 0
    std_gc  = statistics.stdev(kept_gc) if len(kept_gc)>1 else 0

    # summary.tsv
    with open(os.path.join(OUT_DIR, "summary.tsv"), "w") as f:
        f.write("Metric\tValue\n")
        f.write(f"Total\t{total}\n")
        f.write(f"MinGC\t{MIN_GC}\n")
        f.write(f"MaxGC\t{MAX_GC}\n")
        f.write(f"Kept\t{kept}\n")
        f.write(f"FilteredOut\t{total-kept}\n")
        f.write(f"TooLowGC\t{low_gc}\n")
        f.write(f"TooHighGC\t{high_gc}\n")
        f.write(f"MeanGC\t{mean_gc:.3f}\n")
        f.write(f"MedianGC\t{med_gc}\n")
        f.write(f"StdDevGC\t{std_gc:.3f}\n")

    # pretty summary
    with open(os.path.join(OUT_DIR, "summary_pretty.txt"), "w") as f:
        f.write("=============== GC FILTER SUMMARY ===============\n")
        f.write(f"Total introns        : {total}\n")
        f.write(f"Accepted ({MIN_GC}-{MAX_GC}%) : {kept}\n")
        f.write(f"Filtered out         : {total-kept}\n")
        f.write(f"   - GC too low      : {low_gc}\n")
        f.write(f"   - GC too high     : {high_gc}\n")
        f.write("\n")
        f.write(f"Mean GC%             : {mean_gc:.3f}\n")
        f.write(f"Median GC%           : {med_gc}\n")
        f.write(f"Std deviation        : {std_gc:.3f}\n")
        f.write("=================================================\n")

    # ======================================================
    # ORIGINAL PLOTS 
    # ======================================================

    # Histogram
    plt.figure(figsize=(8,5))
    plt.hist(all_gc, bins=40)
    plt.title("GC% Distribution (All Introns)")
    plt.xlabel("GC%")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "histogram.png"))
    plt.close()

    # Boxplot
    plt.figure(figsize=(5,5))
    plt.boxplot(all_gc)
    plt.title("GC% Boxplot")
    plt.ylabel("GC%")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "boxplot.png"))
    plt.close()

    # Barplot
    plt.figure(figsize=(6,5))
    plt.bar(["Low GC","In Range","High GC"], [low_gc, kept, high_gc])
    plt.title("GC Category Distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "barplot.png"))
    plt.close()

    # Scatter
    plt.figure(figsize=(8,5))
    plt.scatter(range(len(all_gc)), all_gc, s=5)
    plt.title("Index vs GC%")
    plt.xlabel("Index")
    plt.ylabel("GC%")
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, "scatter.png"))
    plt.close()

    # ======================================================
    # ENHANCED TERMINAL SUMMARY
    # ======================================================

    print("\n" + "="*65)
    print("GC FILTER SUMMARY - ENHANCED")
    print("="*65)
    print(f"Total introns          : {total:,}")
    print(f"Accepted ({MIN_GC}-{MAX_GC}%)   : {kept:,} ({stats['retention_rate']:.1f}%)")
    print(f"Filtered out           : {total-kept:,} ({stats['filtering_rate']:.1f}%)")
    print(f"   - GC too low        : {low_gc:,} ({stats['low_gc_rate']:.1f}%)")
    print(f"   - GC too high       : {high_gc:,} ({stats['high_gc_rate']:.1f}%)")
    print()
    print("GC STATISTICS (Kept Introns):")
    print(f"  Mean GC%             : {stats['kept_mean']:.2f}%")
    print(f"  Median GC%           : {stats['kept_median']:.2f}%")
    print(f"  Std deviation        : {stats['kept_std']:.2f}%")
    print(f"  Range                : {stats['kept_min']:.2f}% - {stats['kept_max']:.2f}%")
    print()
    print("GC STATISTICS (All Introns):")
    print(f"  Mean GC%             : {stats['all_mean']:.2f}%")
    print(f"  Median GC%           : {stats['all_median']:.2f}%")
    print("="*65)

    print(f"\n[✓] DONE — All GC-filter outputs written to {OUT_DIR}/")
    print(f"[📊] Enhanced GC statistics: {STATS_DIR}/")
    print(f"[📈] Enhanced GC plots: {PLOT_DIR}/")
    print(f"[📋] Original outputs preserved")

# ============================================================
if __name__ == "__main__":
    main()