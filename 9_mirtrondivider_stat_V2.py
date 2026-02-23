#!/usr/bin/env python3
"""
MIRTRON DIVIDER v8 — DROSHA-FREE ONLY + ENHANCED STATISTICS
===========================================================
This script groups ONLY drosha-free = YES introns into:

    • canonical      (O5<=1 and O3<=1)
    • 3p-tailed      (O5<=1 and O3>1)
    • 5p-tailed      (O5>1 and O3<=1)
    • both-tailed    (O5>1 and O3>1)

Inputs:
    drosha_free_results/<source>/<class>_drosha.tsv
    drosha_free_results/<source>/<class>_drosha.fa

Outputs:
    mirtron_divide/<source>/*.tsv
    mirtron_divide/<source>/*.fa
    mirtron_divide/summary_clean.tsv
    mirtron_divide/statistics/*.tsv
    mirtron_divide/statistics/*.txt
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import statistics

print("=== MIRTRON DIVIDER v8 (Drosha-Free ONLY + ENHANCED STATISTICS) ===\n")

# ---------------------------------------------
# INPUT / OUTPUT DIRECTORIES
# ---------------------------------------------
default_in = "drosha_free_results"
in_dir = input(f"➡️  Drosha-free results folder [default: {default_in}]: ").strip() or default_in
BASE_IN  = Path(in_dir)

default_out = "mirtron_divide"
out_dir = input(f"➡️  Output folder [default: {default_out}]: ").strip() or default_out
BASE_OUT = Path(out_dir)

STATS_OUT = BASE_OUT / "statistics"
STATS_OUT.mkdir(exist_ok=True, parents=True)

SOURCES = ["premature", "mature"]
CLASSES = ["very_high", "high", "moderate", "low"]
CATEGORIES = ["canonical", "3p_tailed", "5p_tailed", "both_tailed"]


# ---------------------------------------------
# COLUMN NORMALIZER
# ---------------------------------------------
def normalize_overhang_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    8_overhang_stat_UPDATED.py writes O5_Value / O3_Value.
    Older pipelines may have O5/O3 or Overhang5/Overhang3.
    This function ensures numeric columns named O5 and O3 exist.
    """
    # Prefer explicit value columns if present
    if "O5_Value" in df.columns and "O3_Value" in df.columns:
        df["O5"] = df["O5_Value"]
        df["O3"] = df["O3_Value"]
    else:
        # Try common alternatives
        if "O5" not in df.columns:
            for k in ["Overhang5", "Overhang5_len", "Overhang_5", "Overhang5p", "Overhang5P"]:
                if k in df.columns:
                    df["O5"] = df[k]
                    break
        if "O3" not in df.columns:
            for k in ["Overhang3", "Overhang3_len", "Overhang_3", "Overhang3p", "Overhang3P"]:
                if k in df.columns:
                    df["O3"] = df[k]
                    break
    return df

summary_rows = []
detailed_stats = defaultdict(dict)


def read_fasta(path):
    seqs = {}
    if not path.exists():
        return seqs
    with open(path) as f:
        sid = None
        buff = []
        for line in f:
            if line.startswith(">"):
                if sid:
                    seqs[sid] = "".join(buff)
                sid = line[1:].split("|")[0].strip()
                buff = []
            else:
                buff.append(line.strip())
        if sid:
            seqs[sid] = "".join(buff)
    return seqs


def write_fasta(df, fmap, output):
    with open(output, "w") as out:
        for _, row in df.iterrows():
            cid = row["CoordID"] if "CoordID" in row and pd.notna(row["CoordID"]) else row.get("id", "")
            if cid in fmap:
                out.write(f">{cid}\n{fmap[cid]}\n")


def calculate_category_statistics(df, category_name):
    """Calculate detailed statistics for each mirtron category"""
    if len(df) == 0:
        return {
            'count': 0,
            'o5_mean': 0, 'o5_median': 0, 'o5_std': 0,
            'o3_mean': 0, 'o3_median': 0, 'o3_std': 0,
            'length_mean': 0, 'length_median': 0, 'length_std': 0,
            'gc_mean': 0, 'gc_median': 0
        }
    
    stats = {'count': len(df)}
    
    # Overhang statistics
    o5_values = df["O5"].dropna()
    o3_values = df["O3"].dropna()
    
    if len(o5_values) > 0:
        stats.update({
            'o5_mean': statistics.mean(o5_values),
            'o5_median': statistics.median(o5_values),
            'o5_std': statistics.stdev(o5_values) if len(o5_values) > 1 else 0
        })
    
    if len(o3_values) > 0:
        stats.update({
            'o3_mean': statistics.mean(o3_values),
            'o3_median': statistics.median(o3_values),
            'o3_std': statistics.stdev(o3_values) if len(o3_values) > 1 else 0
        })
    
    # Length statistics
    if "Length" in df.columns:
        length_values = pd.to_numeric(df["Length"], errors='coerce').dropna()
        if len(length_values) > 0:
            stats.update({
                'length_mean': statistics.mean(length_values),
                'length_median': statistics.median(length_values),
                'length_std': statistics.stdev(length_values) if len(length_values) > 1 else 0
            })
    
    # GC content statistics
    if "GC" in df.columns:
        gc_values = pd.to_numeric(df["GC"], errors='coerce').dropna()
        if len(gc_values) > 0:
            stats.update({
                'gc_mean': statistics.mean(gc_values),
                'gc_median': statistics.median(gc_values)
            })
    
    return stats


def analyze_overhang_patterns(df):
    """Analyze detailed overhang patterns and distributions"""
    patterns = {
        'o5_exact_0': len(df[df["O5"] == 0]),
        'o5_exact_1': len(df[df["O5"] == 1]),
        'o5_exact_2': len(df[df["O5"] == 2]),
        'o5_3plus': len(df[df["O5"] >= 3]),
        
        'o3_exact_0': len(df[df["O3"] == 0]),
        'o3_exact_1': len(df[df["O3"] == 1]),
        'o3_exact_2': len(df[df["O3"] == 2]),
        'o3_3plus': len(df[df["O3"] >= 3]),
        
        'both_0': len(df[(df["O5"] == 0) & (df["O3"] == 0)]),
        'perfect_canonical': len(df[(df["O5"] == 1) & (df["O3"] == 1)]),
    }
    
    return patterns


def generate_comprehensive_statistics(detailed_stats):
    """Generate comprehensive statistical reports"""
    
    # Category distribution analysis
    category_report = STATS_OUT / "category_distribution_analysis.tsv"
    with open(category_report, 'w') as f:
        f.write("Source\tClass\tCategory\tCount\tO5_Mean\tO5_Median\tO3_Mean\tO3_Median\tLength_Mean\tGC_Mean\n")
        for (src, cl), categories in detailed_stats.items():
            for category, stats in categories.items():
                if category in CATEGORIES:
                    f.write(f"{src}\t{cl}\t{category}\t{stats['count']}\t"
                           f"{stats.get('o5_mean', 0):.2f}\t{stats.get('o5_median', 0):.1f}\t"
                           f"{stats.get('o3_mean', 0):.2f}\t{stats.get('o3_median', 0):.1f}\t"
                           f"{stats.get('length_mean', 0):.1f}\t{stats.get('gc_mean', 0):.2f}\n")
    
    # Overall summary statistics
    overall_stats = STATS_OUT / "overall_mirtron_statistics.txt"
    with open(overall_stats, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("COMPREHENSIVE MIRTRON CATEGORY STATISTICS\n")
        f.write("=" * 70 + "\n\n")
        
        # Calculate totals
        total_by_category = defaultdict(int)
        total_by_source = defaultdict(int)
        total_by_class = defaultdict(int)
        
        for (src, cl), categories in detailed_stats.items():
            for category, stats in categories.items():
                if category in CATEGORIES:
                    total_by_category[category] += stats['count']
                    total_by_source[src] += stats['count']
                    total_by_class[cl] += stats['count']
        
        # Overall summary
        f.write("OVERALL SUMMARY BY CATEGORY:\n")
        f.write("-" * 50 + "\n")
        total_all = sum(total_by_category.values())
        for category in CATEGORIES:
            count = total_by_category[category]
            percentage = (count / total_all * 100) if total_all > 0 else 0
            f.write(f"{category.replace('_', ' ').title():<15}: {count:>6,} ({percentage:5.1f}%)\n")
        
        f.write(f"\n{'Total':<15}: {total_all:>6,}\n\n")
        
        # Source distribution
        f.write("DISTRIBUTION BY SOURCE:\n")
        f.write("-" * 50 + "\n")
        for source in SOURCES:
            count = total_by_source[source]
            percentage = (count / total_all * 100) if total_all > 0 else 0
            f.write(f"{source.title():<15}: {count:>6,} ({percentage:5.1f}%)\n")
        
        f.write("\n")
        
        # Class distribution
        f.write("DISTRIBUTION BY QUALITY CLASS:\n")
        f.write("-" * 50 + "\n")
        for class_name in CLASSES:
            count = total_by_class[class_name]
            percentage = (count / total_all * 100) if total_all > 0 else 0
            f.write(f"{class_name.replace('_', ' ').title():<15}: {count:>6,} ({percentage:5.1f}%)\n")
    
    # Pattern analysis
    pattern_report = STATS_OUT / "overhang_pattern_analysis.tsv"
    with open(pattern_report, 'w') as f:
        f.write("Source\tClass\tPattern\tCount\n")
        for (src, cl), categories in detailed_stats.items():
            if 'patterns' in categories:
                for pattern, count in categories['patterns'].items():
                    f.write(f"{src}\t{cl}\t{pattern}\t{count}\n")
    
    # Success rate analysis
    success_report = STATS_OUT / "success_rate_analysis.tsv"
    with open(success_report, 'w') as f:
        f.write("Source\tClass\tTotal_Drosha_Free\tCanonical_Pct\t3pTailed_Pct\t5pTailed_Pct\tBothTailed_Pct\n")
        for (src, cl), categories in detailed_stats.items():
            total_drosha_free = sum(stats['count'] for cat, stats in categories.items() if cat in CATEGORIES)
            if total_drosha_free > 0:
                canonical_pct = (categories.get('canonical', {}).get('count', 0) / total_drosha_free * 100)
                tailed_3p_pct = (categories.get('3p_tailed', {}).get('count', 0) / total_drosha_free * 100)
                tailed_5p_pct = (categories.get('5p_tailed', {}).get('count', 0) / total_drosha_free * 100)
                tailed_both_pct = (categories.get('both_tailed', {}).get('count', 0) / total_drosha_free * 100)
                
                f.write(f"{src}\t{cl}\t{total_drosha_free}\t"
                       f"{canonical_pct:.2f}\t{tailed_3p_pct:.2f}\t{tailed_5p_pct:.2f}\t{tailed_both_pct:.2f}\n")


# MAIN LOOP
for src in SOURCES:
    print(f"=== SOURCE: {src.upper()} ===")
    OUT_SRC = BASE_OUT / src
    OUT_SRC.mkdir(exist_ok=True)

    for cl in CLASSES:
        tsv_in = BASE_IN / src / f"{cl}_drosha.tsv"
        fa_in  = BASE_IN / src / f"{cl}_drosha.fa"

        if not tsv_in.exists():
            print(f"[SKIP] {tsv_in} not found.")
            continue

        df = pd.read_csv(tsv_in, sep="\t", dtype=str)

        # Convert O5 + O3 to numeric
        df["O5"] = pd.to_numeric(df["O5"], errors="coerce")
        df["O3"] = pd.to_numeric(df["O3"], errors="coerce")
        
        # Convert other numeric columns if present
        if "Length" in df.columns:
            df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
        if "GC" in df.columns:
            df["GC"] = pd.to_numeric(df["GC"], errors="coerce")

        # ONLY DROSHA-FREE = YES
        df_yes = df[df["Drosha_Free"] == "YES"]
        total_drosha_free = len(df_yes)

        if total_drosha_free == 0:
            print(f"[{src}/{cl}] No Drosha-Free introns found.")
            summary_rows.append({
                "Source": src,
                "Class": cl,
                "Canonical": 0,
                "Tail_3p": 0,
                "Tail_5p": 0,
                "Tail_Both": 0,
                "Total": 0
            })
            continue

        # Load FASTA sequences
        fasta_map = read_fasta(fa_in)

        # CATEGORIES
        canonical     = df_yes[(df_yes["O5"] <= 1) & (df_yes["O3"] <= 1)]
        tailed_3p     = df_yes[(df_yes["O5"] <= 1) & (df_yes["O3"] > 1)]
        tailed_5p     = df_yes[(df_yes["O5"] > 1) & (df_yes["O3"] <= 1)]
        tailed_both   = df_yes[(df_yes["O5"] > 1) & (df_yes["O3"] > 1)]

        # Calculate detailed statistics for each category
        detailed_stats[(src, cl)]['canonical'] = calculate_category_statistics(canonical, 'canonical')
        detailed_stats[(src, cl)]['3p_tailed'] = calculate_category_statistics(tailed_3p, '3p_tailed')
        detailed_stats[(src, cl)]['5p_tailed'] = calculate_category_statistics(tailed_5p, '5p_tailed')
        detailed_stats[(src, cl)]['both_tailed'] = calculate_category_statistics(tailed_both, 'both_tailed')
        
        # Analyze overhang patterns
        detailed_stats[(src, cl)]['patterns'] = analyze_overhang_patterns(df_yes)

        # OUTPUT FILES
        canonical.to_csv(OUT_SRC / f"{cl}_canonical.tsv", sep="\t", index=False)
        tailed_3p.to_csv(OUT_SRC / f"{cl}_3p_tailed.tsv", sep="\t", index=False)
        tailed_5p.to_csv(OUT_SRC / f"{cl}_5p_tailed.tsv", sep="\t", index=False)
        tailed_both.to_csv(OUT_SRC / f"{cl}_both_tailed.tsv", sep="\t", index=False)

        write_fasta(canonical, fasta_map, OUT_SRC / f"{cl}_canonical.fa")
        write_fasta(tailed_3p, fasta_map, OUT_SRC / f"{cl}_3p_tailed.fa")
        write_fasta(tailed_5p, fasta_map, OUT_SRC / f"{cl}_5p_tailed.fa")
        write_fasta(tailed_both, fasta_map, OUT_SRC / f"{cl}_both_tailed.fa")

        summary_rows.append({
            "Source": src,
            "Class": cl,
            "Canonical": len(canonical),
            "Tail_3p": len(tailed_3p),
            "Tail_5p": len(tailed_5p),
            "Tail_Both": len(tailed_both),
            "Total": len(canonical) + len(tailed_3p) + len(tailed_5p) + len(tailed_both)
        })

        print(f"[{src}/{cl}] YES introns → Canonical={len(canonical)}, 3p={len(tailed_3p)}, 5p={len(tailed_5p)}, Both={len(tailed_both)}")


# GENERATE COMPREHENSIVE STATISTICS
print("\n=== GENERATING COMPREHENSIVE STATISTICS ===")
generate_comprehensive_statistics(detailed_stats)

# FINAL SUMMARY TABLE
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(BASE_OUT / "summary_clean.tsv", sep="\t", index=False)

print("\n=== FINAL CLEAN SUMMARY (Drosha-Free ONLY) ===\n")
print(f"{'Source':<12}{'Class':<12}{'Canonical':<12}{'3p':<6}{'5p':<6}{'Both':<6}{'Total':<6}")
print("-" * 70)

for _, r in summary_df.iterrows():
    print(f"{r['Source']:<12}{r['Class']:<12}{int(r['Canonical']):<12}{int(r['Tail_3p']):<6}{int(r['Tail_5p']):<6}{int(r['Tail_Both']):<6}{int(r['Total']):<6}")

print("\n[OK] Saved under mirtron_divide/")
print("📊 Enhanced statistics saved under mirtron_divide/statistics/")
print("=== MIRTRON DIVIDER v8 COMPLETE ===")