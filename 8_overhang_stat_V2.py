#!/usr/bin/env python3
"""
DROSHA-FREE CLASSIFIER (Daneel Ultra v5)
========================================
• Only BLAST matched rows are processed
• Uses O5 / O3 columns
• Realtime logging
• Generates FASTA and TSV for YES hits
• FINAL CLEAN SUMMARY TABLE (publication-ready)
• ENHANCED STATISTICAL ANALYSIS
"""

import csv
import time
import statistics
from pathlib import Path
from datetime import timedelta
from collections import Counter, defaultdict
from Bio import SeqIO

print("=== DROSHA-FREE CLASSIFIER (v5 — BLAST-only, clean summary + STATISTICS) ===")

# ---------------------------------------------
# USER-DEFINED OVERHANG CUTOFF
# ---------------------------------------------
print("Choose maximum allowed Overhang for Drosha-Free detection:")
print("1) 1 nt")
print("2) 2 nt")
print("3) 3 nt")
choice = input("Enter 1, 2, or 3: ").strip()

if choice not in {"1", "2", "3"}:
    raise SystemExit("Invalid choice. Please enter 1, 2, or 3.")

MAX_OVERHANG = int(choice)
print(f"[INFO] Drosha-Free rule: O5 <= {MAX_OVERHANG} and O3 <= {MAX_OVERHANG}\n")

CLASSES = ["very_high", "high", "moderate", "low"]
DB_TYPES = ["premature", "mature"]

# ---------------------------------------------
# INPUT / OUTPUT DIRECTORIES
# ---------------------------------------------
print("📁 Select the BLAST results folder produced by autoblast2.py (or previous BLAST step).")
print("Expected layouts supported:")
print("  A) autoblast2.py: <BASE_IN>/{premature|mature}/{class}.tsv and <BASE_IN>/{premature|mature}/fasta/{class}.fa")
print("  B) legacy:        <BASE_IN>/{premature|mature}/{class}.tsv and <BASE_IN>/{premature|mature}/{class}.fa")
default_in = "blast_results_ultimate"
in_dir = input(f"➡️  BLAST results folder [default: {default_in}]: ").strip() or default_in
BASE_IN = Path(in_dir)

default_out = "drosha_free_results"
out_dir = input(f"➡️  Output folder [default: {default_out}]: ").strip() or default_out
BASE_OUT = Path(out_dir)

# Make output dirs
for db in DB_TYPES:
    (BASE_OUT / db).mkdir(parents=True, exist_ok=True)
(BASE_OUT / "statistics").mkdir(exist_ok=True)


# ---------------------------------------------
# FASTA MAP LOADER
# ---------------------------------------------
def load_fasta(path):
    seqmap = {}
    if not path.exists():
        return seqmap
    for rec in SeqIO.parse(path, "fasta"):
        cid = rec.id.split("|")[0]
        seqmap.setdefault(cid, []).append(str(rec.seq))
    return seqmap


# ---------------------------------------------
# FASTA PATH RESOLVER (supports multiple layouts)
# ---------------------------------------------
def resolve_fasta_path(base_in: Path, db: str, cl: str) -> Path:
    """
    Supports:
      - autoblast2.py layout: base_in/db/fasta/cl.fa
      - legacy layout:        base_in/db/cl.fa
    """
    p1 = base_in / db / "fasta" / f"{cl}.fa"
    if p1.exists():
        return p1
    p2 = base_in / db / f"{cl}.fa"
    if p2.exists():
        return p2
    return p1  # default expected path (even if missing)


# ---------------------------------------------
# STATISTICAL ANALYSIS FUNCTIONS
# ---------------------------------------------
def calculate_detailed_statistics(all_data, max_overhang):
    """Calculate comprehensive statistics for Drosha-Free analysis"""
    stats = {}
    
    # Overall statistics
    total_blast_matches = sum(data['blast_matches'] for data in all_data.values())
    total_drosha_yes = sum(data['drosha_yes'] for data in all_data.values())
    total_drosha_no = sum(data['drosha_no'] for data in all_data.values())
    
    stats['overall'] = {
        'total_blast_matches': total_blast_matches,
        'total_drosha_yes': total_drosha_yes,
        'total_drosha_no': total_drosha_no,
        'overall_yes_percentage': (total_drosha_yes / total_blast_matches * 100) if total_blast_matches > 0 else 0,
        'max_overhang_cutoff': max_overhang
    }
    
    # Database comparison statistics
    db_stats = {}
    for db_type in DB_TYPES:
        db_data = [data for key, data in all_data.items() if key.startswith(db_type)]
        db_blast_matches = sum(data['blast_matches'] for data in db_data)
        db_drosha_yes = sum(data['drosha_yes'] for data in db_data)
        
        db_stats[db_type] = {
            'blast_matches': db_blast_matches,
            'drosha_yes': db_drosha_yes,
            'yes_percentage': (db_drosha_yes / db_blast_matches * 100) if db_blast_matches > 0 else 0,
            'classes_with_hits': len([data for data in db_data if data['blast_matches'] > 0])
        }
    
    stats['database_comparison'] = db_stats
    
    # Class performance statistics
    class_stats = {}
    for class_name in CLASSES:
        class_data = [data for key, data in all_data.items() if key.endswith(class_name)]
        class_blast_matches = sum(data['blast_matches'] for data in class_data)
        class_drosha_yes = sum(data['drosha_yes'] for data in class_data)
        
        class_stats[class_name] = {
            'blast_matches': class_blast_matches,
            'drosha_yes': class_drosha_yes,
            'yes_percentage': (class_drosha_yes / class_blast_matches * 100) if class_blast_matches > 0 else 0,
            'databases_with_hits': len([data for data in class_data if data['blast_matches'] > 0])
        }
    
    stats['class_performance'] = class_stats
    
    # Quality class trend analysis
    quality_trend = []
    for class_name in CLASSES:
        if class_name in class_stats:
            quality_trend.append({
                'class': class_name,
                'yes_percentage': class_stats[class_name]['yes_percentage'],
                'blast_matches': class_stats[class_name]['blast_matches']
            })
    
    stats['quality_trend'] = quality_trend
    
    # Success rate analysis
    high_quality_classes = ['very_high', 'high']
    low_quality_classes = ['moderate', 'low']
    
    high_quality_yes = sum(class_stats[cls]['drosha_yes'] for cls in high_quality_classes if cls in class_stats)
    high_quality_total = sum(class_stats[cls]['blast_matches'] for cls in high_quality_classes if cls in class_stats)
    
    low_quality_yes = sum(class_stats[cls]['drosha_yes'] for cls in low_quality_classes if cls in class_stats)
    low_quality_total = sum(class_stats[cls]['blast_matches'] for cls in low_quality_classes if cls in class_stats)
    
    stats['success_analysis'] = {
        'high_quality_yes_pct': (high_quality_yes / high_quality_total * 100) if high_quality_total > 0 else 0,
        'low_quality_yes_pct': (low_quality_yes / low_quality_total * 100) if low_quality_total > 0 else 0,
        'high_quality_ratio': high_quality_yes / low_quality_yes if low_quality_yes > 0 else float('inf')
    }
    
    return stats

def write_statistical_reports(stats, output_dir):
    """Write comprehensive statistical reports"""
    
    # Main statistics report
    stats_report = output_dir / "statistics" / "detailed_statistics_report.txt"
    with open(stats_report, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("DROSHA-FREE CLASSIFICATION - DETAILED STATISTICAL REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        # Overall Summary
        f.write("OVERALL SUMMARY:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Maximum Overhang Cutoff: {stats['overall']['max_overhang_cutoff']} nt\n")
        f.write(f"Total BLAST Matches Analyzed: {stats['overall']['total_blast_matches']:,}\n")
        f.write(f"Total Drosha-Free Candidates: {stats['overall']['total_drosha_yes']:,}\n")
        f.write(f"Overall Detection Rate: {stats['overall']['overall_yes_percentage']:.2f}%\n\n")
        
        # Database Comparison
        f.write("DATABASE COMPARISON:\n")
        f.write("-" * 50 + "\n")
        for db_type, db_data in stats['database_comparison'].items():
            f.write(f"{db_type.upper()}:\n")
            f.write(f"  BLAST Matches: {db_data['blast_matches']:,}\n")
            f.write(f"  Drosha-Free: {db_data['drosha_yes']:,}\n")
            f.write(f"  Success Rate: {db_data['yes_percentage']:.2f}%\n")
            f.write(f"  Classes with Hits: {db_data['classes_with_hits']}\n\n")
        
        # Class Performance
        f.write("QUALITY CLASS PERFORMANCE:\n")
        f.write("-" * 50 + "\n")
        for class_name, class_data in stats['class_performance'].items():
            f.write(f"{class_name.replace('_', ' ').title()}:\n")
            f.write(f"  BLAST Matches: {class_data['blast_matches']:,}\n")
            f.write(f"  Drosha-Free: {class_data['drosha_yes']:,}\n")
            f.write(f"  Success Rate: {class_data['yes_percentage']:.2f}%\n")
            f.write(f"  Databases with Hits: {class_data['databases_with_hits']}\n\n")
        
        # Success Analysis
        f.write("SUCCESS RATE ANALYSIS:\n")
        f.write("-" * 50 + "\n")
        success_data = stats['success_analysis']
        f.write(f"High Quality Classes Success Rate: {success_data['high_quality_yes_pct']:.2f}%\n")
        f.write(f"Low Quality Classes Success Rate: {success_data['low_quality_yes_pct']:.2f}%\n")
        f.write(f"High/Low Quality Ratio: {success_data['high_quality_ratio']:.2f}\n\n")
        
        # Quality Trend
        f.write("QUALITY TREND ANALYSIS:\n")
        f.write("-" * 50 + "\n")
        for trend in stats['quality_trend']:
            f.write(f"{trend['class'].replace('_', ' ').title():<12}: {trend['yes_percentage']:6.2f}% ({trend['blast_matches']:>6,} matches)\n")
    
    # CSV format for further analysis
    stats_csv = output_dir / "statistics" / "statistical_data.csv"
    with open(stats_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'Subcategory', 'Metric', 'Value'])
        
        # Overall stats
        writer.writerow(['Overall', 'Summary', 'Total_BLAST_Matches', stats['overall']['total_blast_matches']])
        writer.writerow(['Overall', 'Summary', 'Total_Drosha_Free', stats['overall']['total_drosha_yes']])
        writer.writerow(['Overall', 'Summary', 'Overall_Success_Rate', f"{stats['overall']['overall_yes_percentage']:.2f}%"])
        
        # Database stats
        for db_type, db_data in stats['database_comparison'].items():
            writer.writerow(['Database', db_type, 'BLAST_Matches', db_data['blast_matches']])
            writer.writerow(['Database', db_type, 'Drosha_Free', db_data['drosha_yes']])
            writer.writerow(['Database', db_type, 'Success_Rate', f"{db_data['yes_percentage']:.2f}%"])
        
        # Class stats
        for class_name, class_data in stats['class_performance'].items():
            writer.writerow(['Class', class_name, 'BLAST_Matches', class_data['blast_matches']])
            writer.writerow(['Class', class_name, 'Drosha_Free', class_data['drosha_yes']])
            writer.writerow(['Class', class_name, 'Success_Rate', f"{class_data['yes_percentage']:.2f}%"])
    
    print(f"[STATS] Detailed statistical reports saved to {output_dir / 'statistics'}/")



# ---------------------------------------------
# OVERHANG VALUE PARSER (robust column support)
# ---------------------------------------------
def get_overhang_values(row):
    """
    Tries multiple common column names for 5' and 3' overhang.
    Priority:
      1) O5 / O3
      2) Overhang5 / Overhang3
      3) Overhang5_len / Overhang3_len
      4) Overhang_5 / Overhang_3
      5) Overhang5p / Overhang3p
    Returns (o5, o3) as floats; missing => 999.
    """
    keys_5 = ["O5", "Overhang5", "Overhang5_len", "Overhang_5", "Overhang5p", "Overhang5P"]
    keys_3 = ["O3", "Overhang3", "Overhang3_len", "Overhang_3", "Overhang3p", "Overhang3P"]
    def pick(keys):
        for k in keys:
            v = row.get(k, "")
            if v is None:
                continue
            s = str(v).strip()
            if s == "":
                continue
            try:
                return float(s)
            except:
                continue
        return 999.0
    return pick(keys_5), pick(keys_3)

def analyze_overhang_distribution(rows, max_overhang):
    """Analyze overhang distribution patterns (robust overhang parsing)."""
    o5_values = []
    o3_values = []
    both_within_cutoff = 0
    one_within_cutoff = 0
    both_above_cutoff = 0

    for r in rows:
        o5, o3 = get_overhang_values(r)

        # Skip rows with missing overhangs
        if o5 == 999.0 or o3 == 999.0:
            continue

        o5_values.append(o5)
        o3_values.append(o3)

        if o5 <= max_overhang and o3 <= max_overhang:
            both_within_cutoff += 1
        elif o5 <= max_overhang or o3 <= max_overhang:
            one_within_cutoff += 1
        else:
            both_above_cutoff += 1

    distribution = {
        'both_within_cutoff': both_within_cutoff,
        'one_within_cutoff': one_within_cutoff,
        'both_above_cutoff': both_above_cutoff,
        'total_analyzed': len(o5_values),
        'o5_mean': statistics.mean(o5_values) if o5_values else 0,
        'o3_mean': statistics.mean(o3_values) if o3_values else 0,
        'o5_median': statistics.median(o5_values) if o5_values else 0,
        'o3_median': statistics.median(o3_values) if o3_values else 0,
    }

    return distribution


# ---------------------------------------------
# MAIN EXECUTION
# ---------------------------------------------
summary = {}  # For final clean table
detailed_data = {}  # For statistical analysis

for db in DB_TYPES:
    summary[db] = {}

    print(f"\n=== PROCESSING DATABASE: {db.upper()} ===")
    for cl in CLASSES:

        tsv_in = BASE_IN / db / f"{cl}.tsv"
        fa_in = resolve_fasta_path(BASE_IN, db, cl)

        if not tsv_in.exists():
            print(f"[SKIP] {tsv_in} not found")
            continue

        print(f"\n[CLASS] {cl}")

        # Load TSV
        with open(tsv_in) as f:
            rdr = csv.DictReader(f, delimiter="\t")
            all_rows = list(rdr)

        # ------- FILTER: ONLY BLAST MATCHED ROWS -------
        rows = [r for r in all_rows if r.get("BLAST_HIT", "").strip() != ""]
        blast_matches = len(rows)
        print(f"[INFO] BLAST Matches = {blast_matches}")

        if blast_matches == 0:
            summary[db][cl] = (0, 0, 0)
            detailed_data[f"{db}_{cl}"] = {
                'blast_matches': 0,
                'drosha_yes': 0,
                'drosha_no': 0
            }
            print("[WARN] No BLAST hits, skipping")
            continue

        # Analyze overhang distribution
        overhang_dist = analyze_overhang_distribution(rows, MAX_OVERHANG)
        print(f"[STATS] Overhang distribution: Both≤{MAX_OVERHANG}: {overhang_dist['both_within_cutoff']}, "
              f"One≤{MAX_OVERHANG}: {overhang_dist['one_within_cutoff']}, "
              f"Both>{MAX_OVERHANG}: {overhang_dist['both_above_cutoff']}")

        # Load FASTA seqs
        fastamap = load_fasta(fa_in)

        # Process rows
        out_rows = []
        fasta_out = []
        start = time.time()

        for i, r in enumerate(rows, 1):

            cid = r.get("CoordID", r.get("id"))

            # Get overhangs (robust)
            o5, o3 = get_overhang_values(r)

            drosha = (o5 <= MAX_OVERHANG and o3 <= MAX_OVERHANG)

            r2 = r.copy()
            r2["Drosha_Free"] = "YES" if drosha else "NO"
            r2["O5_Value"] = o5
            r2["O3_Value"] = o3
            out_rows.append(r2)

            # FASTA output ONLY if drosha-free
            if drosha and cid in fastamap:
                for seq in fastamap[cid]:
                    fasta_out.append((cid, seq))

            # Realtime
            if i % 25 == 0:
                elapsed = time.time() - start
                rate = i / elapsed if elapsed > 0 else 0
                eta = (blast_matches - i) / rate if rate > 0 else 0
                print(f"[{db}:{cl}] {i}/{blast_matches} | ETA={timedelta(seconds=int(eta))}", end="\r")

        # Output TSV
        tsv_out = BASE_OUT / db / f"{cl}_drosha.tsv"
        with open(tsv_out, "w", newline="") as f:
            fn = list(out_rows[0].keys())
            w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)

        # Output FASTA
        fa_out = BASE_OUT / db / f"{cl}_drosha.fa"
        with open(fa_out, "w") as f:
            for cid, seq in fasta_out:
                f.write(f">{cid}|drosha_free\n{seq}\n")

        # Collect for clean summary
        yes = sum(r["Drosha_Free"] == "YES" for r in out_rows)
        no = sum(r["Drosha_Free"] == "NO" for r in out_rows)

        summary[db][cl] = (blast_matches, yes, no)
        detailed_data[f"{db}_{cl}"] = {
            'blast_matches': blast_matches,
            'drosha_yes': yes,
            'drosha_no': no,
            'overhang_distribution': overhang_dist
        }

# ---------------------------------------------
# STATISTICAL ANALYSIS
# ---------------------------------------------
print("\n\n=== GENERATING COMPREHENSIVE STATISTICS ===")
stats = calculate_detailed_statistics(detailed_data, MAX_OVERHANG)
write_statistical_reports(stats, BASE_OUT)

# ---------------------------------------------
# FINAL CLEAN SUMMARY TABLE (formatted)
# ---------------------------------------------
print("\n\n=== FINAL CLEAN SUMMARY TABLE ===\n")

def pretty(db):
    print(f"### {db.upper()} DATABASE ###\n")
    print(f"{'Class':<12}{'BLAST Matches':<15}{'YES':<10}{'NO':<10}{'YES %':<10}")
    print("-" * 60)

    for cl in CLASSES:
        if cl not in summary[db]:
            continue
        total, yes, no = summary[db][cl]
        pct = (yes / total * 100) if total > 0 else 0
        print(f"{cl:<12}{total:<15}{yes:<10}{no:<10}{pct:>6.1f}%")
    print()

pretty("premature")
pretty("mature")

# Also write to a clean TSV
clean_out = BASE_OUT / "drosha_summary_clean.tsv"
with open(clean_out, "w") as f:
    f.write("DB\tClass\tBLAST_Matches\tDroshaFree_Yes\tDroshaFree_No\tYesPercent\n")
    for db in DB_TYPES:
        for cl in CLASSES:
            if cl not in summary[db]:
                continue
            total, yes, no = summary[db][cl]
            pct = (yes / total * 100) if total > 0 else 0
            f.write(f"{db}\t{cl}\t{total}\t{yes}\t{no}\t{pct:.2f}\n")

print("\n[OK] Summary written to drosha_summary_clean.tsv")
print("\n=== DROSHA-FREE CLASSIFIER v5 COMPLETE ===")
print("📊 Enhanced statistical analysis complete!")
print("📈 Check 'drosha_free_results/statistics/' for detailed reports")