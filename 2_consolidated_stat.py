#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
CONSOLIDATED_ENHANCED (ID SAFE + STATISTICS)
- Adding statistics and analyses
- Detailed reports and graphs
- Input/Output structure stable
"""

import os
import csv
import math
import re
import time
import statistics
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict, Counter

print("=== CONSOLIDATED_ENHANCED (ID SAFE + STATISTICS) ===")
BASE_DIR = os.path.abspath(os.path.dirname(__file__))

INPUT_TSV = os.path.join(BASE_DIR, "intron_outputs", "introns.tsv")
INPUT_FA  = os.path.join(BASE_DIR, "intron_outputs", "introns.fa")
OUT_DIR   = os.path.join(BASE_DIR, "consolidated_introns")

os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# İSTATİSTİK FONKSİYONLARI
# ============================================================
def calculate_comprehensive_stats(records):
    """Kapsamlı istatistikler hesapla"""
    stats = {}
    
    # Basic Statistics
    stats['total_introns'] = len(records)
    stats['total_isoforms'] = sum(r['IsoformCount'] for r in records)
    stats['reduction_ratio'] = stats['total_isoforms'] / stats['total_introns'] if stats['total_introns'] > 0 else 0
    stats['compression_rate'] = ((1 - stats['total_introns']/stats['total_isoforms'])*100) if stats['total_isoforms'] > 0 else 0
    
    # Sequence length statistics
    lengths = [len(r['Sequence']) for r in records if r['Sequence']]
    stats['lengths'] = lengths
    if lengths:
        stats['length_mean'] = statistics.mean(lengths)
        stats['length_median'] = statistics.median(lengths)
        stats['length_std'] = statistics.stdev(lengths) if len(lengths) > 1 else 0
        stats['length_min'] = min(lengths)
        stats['length_max'] = max(lengths)
        stats['length_percentiles'] = {
            p: np.percentile(lengths, p) for p in [25, 50, 75, 90, 95, 99]
        }
    else:
        stats.update({f'length_{k}': 0 for k in ['mean', 'median', 'std', 'min', 'max']})
        stats['length_percentiles'] = {}
    
    # Isoform number statistics
    iso_counts = [r['IsoformCount'] for r in records]
    stats['iso_counts'] = iso_counts
    if iso_counts:
        stats['iso_mean'] = statistics.mean(iso_counts)
        stats['iso_median'] = statistics.median(iso_counts)
        stats['iso_max'] = max(iso_counts)
        stats['iso_percentiles'] = {
            p: np.percentile(iso_counts, p) for p in [50, 75, 90, 95, 99]
        }
    
    # Source distribution
    stats['sources'] = Counter(r['SequenceSource'] for r in records)
    
    return stats

def generate_statistical_plots(records, stats, out_dir):
    """İstatistiksel grafikler oluştur"""
    plot_dir = os.path.join(out_dir, "statistical_plots")
    os.makedirs(plot_dir, exist_ok=True)
    
    # 1. Length Distribution
    plt.figure(figsize=(12, 10))
    
    # 1a. Length Histogram
    plt.subplot(2, 2, 1)
    lengths = stats['lengths']
    if lengths:
        plt.hist(lengths, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        plt.axvline(stats['length_mean'], color='red', linestyle='--', label=f'Mean: {stats["length_mean"]:.1f}')
        plt.axvline(stats['length_median'], color='orange', linestyle='--', label=f'Median: {stats["length_median"]:.1f}')
        plt.xlabel('Intron Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Intron Length Distribution')
        plt.legend()
        plt.yscale('log')
    
    # 1b. Isoform Count Distribution
    plt.subplot(2, 2, 2)
    iso_counts = stats['iso_counts']
    if iso_counts:
        plt.hist(iso_counts, bins=min(20, max(iso_counts)), alpha=0.7, color='lightgreen', edgecolor='black')
        plt.xlabel('Isoforms per Intron')
        plt.ylabel('Frequency')
        plt.title('Isoform Count Distribution')
        plt.yscale('log')
    
    # 1c. Source Distribution
    plt.subplot(2, 2, 3)
    source_data = stats['sources']
    if source_data:
        labels, values = zip(*source_data.items())
        plt.bar(labels, values, color='orange', alpha=0.7, edgecolor='black')
        plt.xlabel('Sequence Source')
        plt.ylabel('Count')
        plt.title('Sequence Source Distribution')
        plt.xticks(rotation=45)
    
    # 1d. Length vs Isoform Count
    plt.subplot(2, 2, 4)
    if lengths and iso_counts:
        plt.scatter(lengths, iso_counts, alpha=0.6, s=10)
        plt.xlabel('Intron Length (bp)')
        plt.ylabel('Isoform Count')
        plt.title('Intron Length vs Isoform Count')
        plt.yscale('log')
        plt.xscale('log')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'comprehensive_statistics.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"[STATS] Statistical plots saved to {plot_dir}/")

def save_detailed_statistics_report(stats, records, out_dir):
    """Detaylı istatistik raporu oluştur"""
    report_path = os.path.join(out_dir, "detailed_statistics_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("CONSOLIDATED INTRONS - DETAILED STATISTICAL REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        # General summary
        f.write("OVERVIEW SUMMARY:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total Introns Consolidated: {stats['total_introns']:,}\n")
        f.write(f"Total Isoforms Processed: {stats['total_isoforms']:,}\n")
        f.write(f"Reduction Ratio: {stats['reduction_ratio']:.2f} isoforms/intron\n")
        f.write(f"Compression Rate: {stats['compression_rate']:.1f}%\n")
        f.write(f"Average Isoforms per Intron: {stats['iso_mean']:.2f}\n\n")
        
        # Length Statistics
        f.write("LENGTH STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats['lengths']:
            f.write(f"Mean Length: {stats['length_mean']:.2f} bp\n")
            f.write(f"Median Length: {stats['length_median']:.2f} bp\n")
            f.write(f"Standard Deviation: {stats['length_std']:.2f} bp\n")
            f.write(f"Range: {stats['length_min']:,} - {stats['length_max']:,} bp\n\n")
            
            f.write("LENGTH PERCENTILES:\n")
            for p, value in stats['length_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.2f} bp\n")
            f.write("\n")
            
            # Length categories
            len_cats = {
                "Very Short (<50 bp)": len([l for l in stats['lengths'] if l < 50]),
                "Short (50-100 bp)": len([l for l in stats['lengths'] if 50 <= l < 100]),
                "Medium (100-500 bp)": len([l for l in stats['lengths'] if 100 <= l < 500]),
                "Long (500-2000 bp)": len([l for l in stats['lengths'] if 500 <= l < 2000]),
                "Very Long (≥2000 bp)": len([l for l in stats['lengths'] if l >= 2000])
            }
            
            f.write("LENGTH CATEGORIES:\n")
            for cat, count in len_cats.items():
                if count > 0:
                    pct = (count / len(stats['lengths'])) * 100
                    f.write(f"  {cat}: {count:,} ({pct:.1f}%)\n")
            f.write("\n")
        
        # Isoform Statistics
        f.write("ISOFORM COUNT STATISTICS:\n")
        f.write("-" * 50 + "\n")
        if stats['iso_counts']:
            f.write(f"Mean Isoforms per Intron: {stats['iso_mean']:.2f}\n")
            f.write(f"Median Isoforms per Intron: {stats['iso_median']:.2f}\n")
            f.write(f"Maximum Isoforms: {stats['iso_max']}\n\n")
            
            f.write("ISOFORM PERCENTILES:\n")
            for p, value in stats['iso_percentiles'].items():
                f.write(f"  {p}th percentile: {value:.2f}\n")
            f.write("\n")
            
            # Isoform categories
            iso_cats = {
                "Single Isoform": len([c for c in stats['iso_counts'] if c == 1]),
                "2-5 Isoforms": len([c for c in stats['iso_counts'] if 2 <= c <= 5]),
                "6-10 Isoforms": len([c for c in stats['iso_counts'] if 6 <= c <= 10]),
                "11-20 Isoforms": len([c for c in stats['iso_counts'] if 11 <= c <= 20]),
                ">20 Isoforms": len([c for c in stats['iso_counts'] if c > 20])
            }
            
            f.write("ISOFORM CATEGORIES:\n")
            for cat, count in iso_cats.items():
                if count > 0:
                    pct = (count / len(stats['iso_counts'])) * 100
                    f.write(f"  {cat}: {count:,} ({pct:.1f}%)\n")
            f.write("\n")
        
        # Source Distribution
        f.write("SEQUENCE SOURCE DISTRIBUTION:\n")
        f.write("-" * 50 + "\n")
        for source, count in stats['sources'].items():
            pct = (count / len(records)) * 100
            f.write(f"  {source}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
    
    print(f"[STATS] Detailed statistics report saved to {report_path}")

# ============================================================
# ID SAFETY
# ============================================================
def extract_clean_coord(raw):
    """
    Extract clean CoordID:
    Examples:
    'NC_081805.1:247209-247321:+|rna-XM...'  → 'NC_081805.1:247209-247321:+'
    'chr1:100-200:+'                         → 'chr1:100-200:+'
    """
    pattern = r'([^:|]+:[0-9]+-[0-9]+:[+-])'
    m = re.search(pattern, raw)
    return m.group(1) if m else raw

def make_fasta_header(coord, iso_count, seq_len, src):
    return f">{coord}|isoforms={iso_count}|len={seq_len}|src={src}"

# ============================================================
# CONFIGURATION
# ============================================================
LOG_FREQUENCY = 2000  # ✅ Configurable: Her kaç intronda bir log yazılsın

# ============================================================
# LOAD introns.tsv (ID SAFE)
# ============================================================
if not os.path.exists(INPUT_TSV):
    raise FileNotFoundError(f"[ERROR] Missing: {INPUT_TSV}")

rows = []
print("[INFO] Loading", INPUT_TSV)
with open(INPUT_TSV) as f:
    rdr = csv.DictReader(f, delimiter="\t")
    for r in rdr:
        if "ID" not in r:
            raise ValueError("[ERROR] TSV missing 'ID' column!")
        if "sequence" not in r:
            raise ValueError("[ERROR] TSV missing 'sequence' column!")

        # ✅ ID GÜVENLİĞİ: CoordID'yi temizle
        clean = extract_clean_coord(r["ID"])
        r["CoordID"] = clean
        r["Sequence"] = r["sequence"]
        rows.append(r)

if not rows:
    raise ValueError("[ERROR] introns.tsv empty")

print(f"[INFO] Loaded {len(rows)} introns")

# ============================================================
# LOAD introns.fa (ID SAFE)
# ============================================================
seq_dict = {}
if os.path.exists(INPUT_FA):
    print("[INFO] Loading", INPUT_FA)
    cur = None
    buf = []
    with open(INPUT_FA) as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    # ✅ ID GÜVENLİĞİ: FASTA ID'yi temizle
                    coord = extract_clean_coord(cur)
                    seq_dict[coord] = "".join(buf)
                cur = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if cur:
            # ✅ ID GÜVENLİĞİ: Son entry'i de temizle
            coord = extract_clean_coord(cur)
            seq_dict[coord] = "".join(buf)

else:
    print("[WARN] FASTA missing:", INPUT_FA)

# ============================================================
# CONSOLIDATE — WITH REALTIME LOGGING + ETA
# ============================================================
group = defaultdict(list)

print(f"\n[STEP] Consolidating introns (realtime - log every {LOG_FREQUENCY})...\n")

start_time = time.time()  # ✅ ETA için başlangıç zamanı

for idx, r in enumerate(rows, 1):
    coord = r["CoordID"]
    rid   = r["ID"]
    seq_tsv = r["Sequence"].strip()

    # ✅ REALTIME LOG WITH ETA
    if idx % LOG_FREQUENCY == 1:   
        elapsed = time.time() - start_time
        progress_pct = (idx / len(rows)) * 100
        
        # ETA calculation
        if idx > 1 and elapsed > 0:
            items_per_second = idx / elapsed
            remaining_items = len(rows) - idx
            eta_seconds = remaining_items / items_per_second if items_per_second > 0 else 0
            
            # Format ETA
            if eta_seconds < 60:
                eta_str = f"{eta_seconds:.0f}s"
            elif eta_seconds < 3600:
                eta_str = f"{eta_seconds/60:.1f}m"
            else:
                eta_str = f"{eta_seconds/3600:.1f}h"
                
            print(f"  → {idx:,}/{len(rows):,} ({progress_pct:.1f}%) - ETA: {eta_str} - {coord}")
        else:
            print(f"  → {idx:,}/{len(rows):,} ({progress_pct:.1f}%) - {coord}")

    # ✅ ID SAFETY: ALL LOOKUPS CLEAN WITH CoordID 
    if seq_tsv:
        group[coord].append((rid, seq_tsv))
        continue

    if coord in seq_dict:
        group[coord].append((rid, seq_dict[coord]))
        continue

    # PREFIX MATCH (clean CoordID)
    hit = ""
    for k, v in seq_dict.items():
        if k.startswith(coord):
            hit = v
            break

    group[coord].append((rid, hit))

consolidation_time = time.time() - start_time
print(f"[INFO] Consolidation grouping complete in {consolidation_time:.1f}s")

# ============================================================
# FINALIZE RECORDS — WITH ETA
# ============================================================
records = []
print(f"\n[STEP] Selecting best sequence per intron (realtime - log every {LOG_FREQUENCY})...\n")

finalize_start = time.time()

for idx, (coord, items) in enumerate(sorted(group.items()), 1):
    # ✅ REALTIME LOG FOR FINALIZATION
    if idx % LOG_FREQUENCY == 1:
        elapsed_final = time.time() - finalize_start
        progress_pct = (idx / len(group)) * 100
        
        if idx > 1 and elapsed_final > 0:
            items_per_second = idx / elapsed_final
            remaining_items = len(group) - idx
            eta_seconds = remaining_items / items_per_second if items_per_second > 0 else 0
            
            if eta_seconds < 60:
                eta_str = f"{eta_seconds:.0f}s"
            elif eta_seconds < 3600:
                eta_str = f"{eta_seconds/60:.1f}m"
            else:
                eta_str = f"{eta_seconds/3600:.1f}h"
                
            print(f"  → {idx:,}/{len(group):,} ({progress_pct:.1f}%) - ETA: {eta_str} - {coord}")
        else:
            print(f"  → {idx:,}/{len(group):,} ({progress_pct:.1f}%) - {coord}")

    non_empty = [s for (_, s) in items if s]
    best_seq = max(non_empty, key=len) if non_empty else ""

    iso_ids = [rid for (rid, _) in items]

    seq_source = "NONE" if not best_seq else "TSV_OR_FASTA"

    records.append({
        "CoordID": coord,
        "Sequence": best_seq,
        "SequenceLength": (len(best_seq) if best_seq else ""),
        "IsoformCount": len(iso_ids),
        "IsoformList": ";".join(iso_ids),
        "SequenceSource": seq_source
    })

finalize_time = time.time() - finalize_start
print(f"[INFO] Finalization complete in {finalize_time:.1f}s")

# ============================================================
# WRITE OUTPUTS (ID SAFE FASTA HEADERS)
# ============================================================
print("\n[STEP] Writing consolidated outputs...")

out_tsv = os.path.join(OUT_DIR, "consolidated.tsv")
with open(out_tsv, "w", newline="") as f:
    fn = ["CoordID","Sequence","SequenceLength",
          "IsoformCount","IsoformList","SequenceSource"]
    w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
    w.writeheader()
    w.writerows(records)

out_fa = os.path.join(OUT_DIR, "consolidated.fa")
with open(out_fa, "w") as f:
    for r in records:
        if r["Sequence"]:
            # ✅ ID GÜVENLİĞİ: Standardize FASTA header
            header = make_fasta_header(
                r["CoordID"], r["IsoformCount"],
                r["SequenceLength"], r["SequenceSource"]
            )
            f.write(header + "\n")
            seq = r["Sequence"]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

print("[INFO] CORE consolidation complete.")

# ============================================================
# PIPELINE SAFETY VALIDATION
# ============================================================
print("\n" + "="*60)
print("PIPELINE SAFETY VALIDATION")
print("="*60)

tsv_coords = {r["CoordID"] for r in records if r["Sequence"]}
fasta_coords = set()

with open(out_fa) as f:
    for line in f:
        if line.startswith(">"):
            fasta_coords.add(line.split("|")[0][1:])

mismatch = tsv_coords - fasta_coords

print(f"TSV CoordIDs:  {len(tsv_coords)}")
print(f"FASTA CoordIDs:{len(fasta_coords)}")
print(f"Perfect match: {len(tsv_coords & fasta_coords)}")
print(f"MISMATCH:      {len(mismatch)}")

if mismatch:
    print("🚨 CRITICAL: FASTA and TSV IDs do NOT match!")
else:
    print("✅ SAFE: FASTA and TSV IDs are perfectly synchronized.")
print("="*60)

# ============================================================
# ENHANCED STATISTICS GENERATION
# ============================================================
print("\n" + "="*60)
print("GENERATING ENHANCED STATISTICS")
print("="*60)

# Calculate comprehensive statistics
stats = calculate_comprehensive_stats(records)

# Generate statistical plots
generate_statistical_plots(records, stats, OUT_DIR)

# Saved detailed report
save_detailed_statistics_report(stats, records, OUT_DIR)

# ============================================================
# SUMMARY STATISTICS 
# ============================================================
print("\n" + "="*60)
print("CONSOLIDATION SUMMARY - ENHANCED")
print("="*60)

print(f"Total isoforms processed: {stats['total_isoforms']:,}")
print(f"Consolidated introns: {stats['total_introns']:,}")
print(f"Reduction ratio: {stats['reduction_ratio']:.2f} isoforms per intron")
print(f"Compression: {stats['compression_rate']:.1f}% reduction")
print(f"Average isoforms per intron: {stats['iso_mean']:.2f}")

# Sequence length summary
if stats['lengths']:
    print(f"\n📊 LENGTH STATISTICS:")
    print(f"  Mean: {stats['length_mean']:.1f} bp")
    print(f"  Median: {stats['length_median']:.1f} bp")
    print(f"  Range: {stats['length_min']:,} - {stats['length_max']:,} bp")
    print(f"  Std Dev: {stats['length_std']:.1f} bp")

# Source summary
print(f"\n📁 SEQUENCE SOURCE DISTRIBUTION:")
for source, count in stats['sources'].items():
    percentage = (count / len(records)) * 100
    print(f"  {source:20}: {count:6,} ({percentage:5.1f}%)")

print("="*60)

total_time = time.time() - start_time
print(f"\n=== CONSOLIDATED v11.1_DANEEL_ENHANCED COMPLETE in {total_time:.1f}s ===")
print(f"📁 Outputs: {OUT_DIR}/")
print(f"📊 Statistics: {OUT_DIR}/detailed_statistics_report.txt")
print(f"📈 Plots: {OUT_DIR}/statistical_plots/")
print("ℹ️  NOTE: PLUGIN A and B (Genome Input & Heatmap) have been removed")