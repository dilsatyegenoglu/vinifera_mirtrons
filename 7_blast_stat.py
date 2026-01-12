#!/usr/bin/env python3
"""
BLAST_DANEEL_ULTRA v2
=====================
FASTA → BLAST → FULL TSV merge → FASTA matched → realtime → ENHANCED STATISTICS

Inputs:
    validator_results/classes/*.fa
    validator_results/classes/*.tsv
    premature.fasta
    mature.fasta

Outputs:
    blast_results/premature/*.tsv
    blast_results/premature/*.fa
    blast_results/mature/*.tsv
    blast_results/mature/*.fa
    blast_results/premature/statistics/*
    blast_results/mature/statistics/*
    blast_results/comparison/*
    summary.txt
    summary.tsv
"""

import csv
import time
import subprocess
import statistics
from pathlib import Path
from collections import Counter, defaultdict
from Bio import SeqIO
from datetime import timedelta

print("=== BLAST_DANEEL_ULTRA v2 (FASTA → TSV merge → FASTA + ENHANCED STATISTICS) ===")

# -----------------------------------
# USER PARAMETERS
# -----------------------------------
pident_min = float(input("Minimum percent identity (e.g. 90): ").strip())
e_min = float(input("Minimum E-value MIN (e.g. 1e-10): ").strip())
e_max = float(input("Maximum E-value MAX (e.g. 1e-2): ").strip())

CATS = ["very_high", "high", "moderate", "low"]
DB_TYPES = ["premature", "mature"]

CLASS_DIR = Path("validator_results/classes")
OUT_DIR = Path("blast_results")
DB_DIR = OUT_DIR / "db"

(DB_DIR).mkdir(parents=True, exist_ok=True)
(OUT_DIR / "premature").mkdir(parents=True, exist_ok=True)
(OUT_DIR / "mature").mkdir(parents=True, exist_ok=True)

HEADER_BLAST = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]

# -----------------------------------
# STATISTICS FUNCTIONS
# -----------------------------------
def calculate_blast_statistics(all_hits, db_type, total_queries):
    """Calculate comprehensive BLAST statistics"""
    stats = {}
    
    if not all_hits:
        return create_empty_stats(db_type, total_queries)
    
    # Extract numerical values
    pidents = [float(hit[2]) for hit in all_hits]
    evalues = [float(hit[10]) for hit in all_hits]
    bitscores = [float(hit[11]) for hit in all_hits]
    lengths = [int(hit[3]) for hit in all_hits]
    
    # General statistics
    stats['general'] = {
        'total_queries': total_queries,
        'total_hits': len(all_hits),
        'queries_with_hits': len(set(hit[0] for hit in all_hits)),
        'hit_percentage': (len(set(hit[0] for hit in all_hits)) / total_queries) * 100 if total_queries > 0 else 0,
        'average_hits_per_query': len(all_hits) / total_queries if total_queries > 0 else 0,
        'max_hits_per_query': max(Counter(hit[0] for hit in all_hits).values()) if all_hits else 0
    }
    
    # Quality metrics
    stats['quality'] = {
        'mean_pident': statistics.mean(pidents) if pidents else 0,
        'median_pident': statistics.median(pidents) if pidents else 0,
        'mean_evalue': statistics.mean(evalues) if evalues else 0,
        'median_evalue': statistics.median(evalues) if evalues else 0,
        'mean_bitscore': statistics.mean(bitscores) if bitscores else 0,
        'median_bitscore': statistics.median(bitscores) if bitscores else 0,
        'mean_length': statistics.mean(lengths) if lengths else 0,
        'median_length': statistics.median(lengths) if lengths else 0
    }
    
    # E-value distribution
    evalue_bins = {
        '1e-50_less': len([e for e in evalues if e < 1e-50]),
        '1e-50_to_1e-20': len([e for e in evalues if 1e-50 <= e < 1e-20]),
        '1e-20_to_1e-10': len([e for e in evalues if 1e-20 <= e < 1e-10]),
        '1e-10_to_1e-5': len([e for e in evalues if 1e-10 <= e < 1e-5]),
        '1e-5_to_1e-2': len([e for e in evalues if 1e-5 <= e <= 1e-2])
    }
    stats['evalue_distribution'] = evalue_bins
    
    # Percent identity distribution
    pident_bins = {
        '95_100': len([p for p in pidents if 95 <= p <= 100]),
        '90_95': len([p for p in pidents if 90 <= p < 95]),
        '85_90': len([p for p in pidents if 85 <= p < 90]),
        '80_85': len([p for p in pidents if 80 <= p < 85]),
        'less_80': len([p for p in pidents if p < 80])
    }
    stats['pident_distribution'] = pident_bins
    
    # Bitscore distribution
    if bitscores:
        bs_min, bs_max = min(bitscores), max(bitscores)
        bs_range = bs_max - bs_min
        bs_bins = {}
        for i in range(5):
            lower = bs_min + (i * bs_range / 5)
            upper = bs_min + ((i + 1) * bs_range / 5)
            bs_bins[f'{lower:.1f}_{upper:.1f}'] = len([b for b in bitscores if lower <= b < upper])
        stats['bitscore_distribution'] = bs_bins
    
    # Multiple hits analysis
    query_hit_count = Counter(hit[0] for hit in all_hits)
    stats['multiple_hits'] = {
        'single_hit_queries': len([q for q, count in query_hit_count.items() if count == 1]),
        'multiple_hit_queries': len([q for q, count in query_hit_count.items() if count > 1]),
        'max_hits_single_query': max(query_hit_count.values()) if query_hit_count else 0
    }
    
    # Top hits (best bitscore per query)
    top_hits = {}
    for hit in all_hits:
        qseqid = hit[0]
        bitscore = float(hit[11])
        if qseqid not in top_hits or bitscore > float(top_hits[qseqid][11]):
            top_hits[qseqid] = hit
    
    stats['top_hits'] = {
        'count': len(top_hits),
        'mean_pident': statistics.mean([float(hit[2]) for hit in top_hits.values()]) if top_hits else 0,
        'mean_evalue': statistics.mean([float(hit[10]) for hit in top_hits.values()]) if top_hits else 0,
        'mean_bitscore': statistics.mean([float(hit[11]) for hit in top_hits.values()]) if top_hits else 0
    }
    
    return stats

def create_empty_stats(db_type, total_queries):
    """Create empty statistics when no hits are found"""
    return {
        'general': {
            'total_queries': total_queries,
            'total_hits': 0,
            'queries_with_hits': 0,
            'hit_percentage': 0,
            'average_hits_per_query': 0,
            'max_hits_per_query': 0
        },
        'quality': {
            'mean_pident': 0, 'median_pident': 0,
            'mean_evalue': 0, 'median_evalue': 0,
            'mean_bitscore': 0, 'median_bitscore': 0,
            'mean_length': 0, 'median_length': 0
        },
        'evalue_distribution': {k: 0 for k in ['1e-50_less', '1e-50_to_1e-20', '1e-20_to_1e-10', '1e-10_to_1e-5', '1e-5_to_1e-2']},
        'pident_distribution': {k: 0 for k in ['95_100', '90_95', '85_90', '80_85', 'less_80']},
        'multiple_hits': {'single_hit_queries': 0, 'multiple_hit_queries': 0, 'max_hits_single_query': 0},
        'top_hits': {'count': 0, 'mean_pident': 0, 'mean_evalue': 0, 'mean_bitscore': 0}
    }

def write_blast_statistics(stats, output_dir, db_type, category_stats=None):
    """Write comprehensive BLAST statistics to files"""
    stats_dir = output_dir / "statistics"
    stats_dir.mkdir(exist_ok=True)
    
    # General summary
    with open(stats_dir / "blast_summary.tsv", "w") as f:
        f.write("Metric\tValue\n")
        for metric, value in stats['general'].items():
            f.write(f"{metric}\t{value}\n")
    
    # Quality metrics
    with open(stats_dir / "quality_metrics.tsv", "w") as f:
        f.write("Metric\tValue\n")
        for metric, value in stats['quality'].items():
            f.write(f"{metric}\t{value}\n")
    
    # E-value distribution
    with open(stats_dir / "evalue_distribution.tsv", "w") as f:
        f.write("Evalue_Range\tCount\n")
        for evalue_range, count in stats['evalue_distribution'].items():
            f.write(f"{evalue_range}\t{count}\n")
    
    # Percent identity distribution
    with open(stats_dir / "pident_distribution.tsv", "w") as f:
        f.write("PIdent_Range\tCount\n")
        for pident_range, count in stats['pident_distribution'].items():
            f.write(f"{pident_range}\t{count}\n")
    
    # Multiple hits analysis
    with open(stats_dir / "multiple_hits.tsv", "w") as f:
        f.write("Metric\tValue\n")
        for metric, value in stats['multiple_hits'].items():
            f.write(f"{metric}\t{value}\n")
    
    # Top hits summary
    with open(stats_dir / "top_hits_summary.tsv", "w") as f:
        f.write("Metric\tValue\n")
        for metric, value in stats['top_hits'].items():
            f.write(f"{metric}\t{value}\n")
    
    # Category statistics if available
    if category_stats:
        with open(stats_dir / "category_stats.tsv", "w") as f:
            f.write("Category\tTotalQueries\tHits\tHitPercentage\tMeanPident\tMeanEvalue\tMeanBitscore\n")
            for cat, cat_data in category_stats.items():
                f.write(f"{cat}\t{cat_data['total_queries']}\t{cat_data['hits']}\t{cat_data['hit_percentage']:.2f}\t{cat_data['mean_pident']:.2f}\t{cat_data['mean_evalue']:.2e}\t{cat_data['mean_bitscore']:.2f}\n")

def create_comparison_statistics(premature_stats, mature_stats, output_dir):
    """Create comparison statistics between premature and mature databases"""
    comp_dir = output_dir / "comparison"
    comp_dir.mkdir(exist_ok=True)
    
    # Database comparison
    with open(comp_dir / "database_comparison.tsv", "w") as f:
        f.write("Metric\tPremature\tMature\tDifference\n")
        
        # General metrics comparison
        general_metrics = ['total_hits', 'queries_with_hits', 'hit_percentage', 'average_hits_per_query']
        for metric in general_metrics:
            pre_val = premature_stats['general'][metric]
            mat_val = mature_stats['general'][metric]
            diff = pre_val - mat_val
            f.write(f"{metric}\t{pre_val}\t{mat_val}\t{diff}\n")
        
        # Quality metrics comparison
        quality_metrics = ['mean_pident', 'mean_evalue', 'mean_bitscore']
        for metric in quality_metrics:
            pre_val = premature_stats['quality'][metric]
            mat_val = mature_stats['quality'][metric]
            diff = pre_val - mat_val
            f.write(f"{metric}\t{pre_val}\t{mat_val}\t{diff}\n")

# -----------------------------------
# BUILD DB + CHECK
# -----------------------------------
def build_db(fasta, dbname):
    subprocess.run([
        "makeblastdb",
        "-in", fasta,
        "-dbtype", "nucl",
        "-out", str(dbname)
    ], check=True)

    # kontrol
    if not Path(str(dbname) + ".nsq").exists():
        raise SystemExit(f"[ERROR] Database creation failed: {dbname}")

print("\n[DB] Building databases...")
build_db("premature.fasta", DB_DIR / "db_premature")
build_db("mature.fasta",    DB_DIR / "db_mature")
print("[OK] Databases ready.\n")

# -----------------------------------
# BLAST RUNNER
# -----------------------------------
def run_blast(query_fa, dbpath):
    cmd = [
        "blastn", "-task", "blastn-short",
        "-query", str(query_fa),
        "-db", str(dbpath),
        "-outfmt", "6 " + " ".join(HEADER_BLAST),
        "-perc_identity", str(pident_min),
        "-evalue", str(e_max)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    lines = result.stdout.strip().split("\n")
    if lines == [''] or len(lines) == 0:
        return []
    return lines

# -----------------------------------
# FASTA PARSER (safe)
# -----------------------------------
def parse_fasta_safe(fapath):
    seqmap = {}  # CoordID → list of sequences
    for rec in SeqIO.parse(fapath, "fasta"):
        header = rec.id
        coord = header.split("|")[0]  # ONLY CoordID
        seqmap.setdefault(coord, []).append(str(rec.seq))
    return seqmap

# -----------------------------------
# MAIN LOOP
# -----------------------------------
summary = []
all_hits_by_db = {db_type: [] for db_type in DB_TYPES}
category_stats_by_db = {db_type: defaultdict(dict) for db_type in DB_TYPES}

for cat in CATS:
    fasta = CLASS_DIR / f"{cat}.fa"
    tsv = CLASS_DIR / f"{cat}.tsv"

    if not fasta.exists() or not tsv.exists():
        print(f"[SKIP] Missing {cat} (FA/TSV)")
        continue

    print(f"\n=== CLASS: {cat} ===")
    seqs = parse_fasta_safe(fasta)

    # Load validator TSV rows
    with open(tsv) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        rows = list(rdr)

    for dbtype in DB_TYPES:
        print(f"\n[{cat} → {dbtype}] Running BLAST...")

        dbpath = DB_DIR / f"db_{dbtype}"
        out_tsv = OUT_DIR / dbtype / f"{cat}.tsv"
        out_fa  = OUT_DIR / dbtype / f"{cat}.fa"

        # Create temp query
        tempfa = Path("temp_blast.fa")
        with open(tempfa, "w") as tf:
            for r in rows:
                cid = r["CoordID"]
                if cid in seqs:
                    for i, s in enumerate(seqs[cid]):
                        tf.write(f">{cid}\n{s}\n")

        # Run BLAST
        raw = run_blast(tempfa, dbpath)

        # Filter
        filtered = []
        for line in raw:
            parts = line.split("\t")
            if len(parts) < 12:
                continue
            evalue = float(parts[10])
            if e_min <= evalue <= e_max:
                filtered.append(parts)

        # Store hits for statistics
        all_hits_by_db[dbtype].extend(filtered)
        
        # Calculate category-specific statistics
        if filtered:
            pidents = [float(hit[2]) for hit in filtered]
            evalues = [float(hit[10]) for hit in filtered]
            bitscores = [float(hit[11]) for hit in filtered]
            
            category_stats_by_db[dbtype][cat] = {
                'total_queries': len(rows),
                'hits': len(filtered),
                'hit_percentage': (len(set(hit[0] for hit in filtered)) / len(rows)) * 100 if rows else 0,
                'mean_pident': statistics.mean(pidents) if pidents else 0,
                'mean_evalue': statistics.mean(evalues) if evalues else 0,
                'mean_bitscore': statistics.mean(bitscores) if bitscores else 0
            }
        else:
            category_stats_by_db[dbtype][cat] = {
                'total_queries': len(rows),
                'hits': 0,
                'hit_percentage': 0,
                'mean_pident': 0,
                'mean_evalue': 0,
                'mean_bitscore': 0
            }

        # REALTIME LOGGING
        total = len(rows)
        start = time.time()
        for i, r in enumerate(rows, 1):
            if i % 50 == 0:
                elapsed = time.time() - start
                rate = i / elapsed if elapsed > 0 else 0
                eta = (total - i) / rate if rate > 0 else 0
                print(f"[{cat}→{dbtype}] {i}/{total} "
                      f"| matched={len(filtered)} "
                      f"| ETA {str(timedelta(seconds=int(eta)))}", end="\r")

        # WRITE TSV (FULL MERGE)
        print(f"\n[WRITE] {out_tsv}")
        with open(out_tsv, "w", newline="") as f:
            fn = list(rows[0].keys()) + [
                "BLAST_HIT","pident","evalue","bitscore",
                "qstart","qend","sstart","send","DB"
            ]
            w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
            w.writeheader()

            hitmap = {}  # CoordID → list of BLAST rows
            for h in filtered:
                cid = h[0]
                hitmap.setdefault(cid, []).append(h)

            for r in rows:
                cid = r["CoordID"]
                if cid in hitmap:
                    for h in hitmap[cid]:
                        _, sid, pid, ln, mis, gap, qs, qe, ss, se, ev, bit = h
                        new = r.copy()
                        new.update({
                            "BLAST_HIT": sid,
                            "pident": pid,
                            "evalue": ev,
                            "bitscore": bit,
                            "qstart": qs,
                            "qend": qe,
                            "sstart": ss,
                            "send": se,
                            "DB": dbtype
                        })
                        w.writerow(new)
                else:
                    new = r.copy()
                    new.update({
                        "BLAST_HIT": "",
                        "pident": "",
                        "evalue": "",
                        "bitscore": "",
                        "qstart": "",
                        "qend": "",
                        "sstart": "",
                        "send": "",
                        "DB": dbtype
                    })
                    w.writerow(new)

        # WRITE FASTA (matched only)
        print(f"[WRITE] {out_fa}")
        with open(out_fa, "w") as f:
            for h in filtered:
                cid = h[0]
                sid = h[1]
                pid = h[2]
                ev = h[10]
                bit = h[11]

                if cid in seqs:
                    for seq in seqs[cid]:
                        f.write(f">{cid}|match={sid}|pid={pid}|ev={ev}|bit={bit}\n")
                        f.write(seq + "\n")

        summary.append((cat, dbtype, len(filtered)))

# -----------------------------------
# COMPREHENSIVE STATISTICS
# -----------------------------------
print("\n=== GENERATING COMPREHENSIVE STATISTICS ===")

# Calculate total queries per database
total_queries_by_db = {}
for dbtype in DB_TYPES:
    total_queries = 0
    for cat in CATS:
        tsv_path = CLASS_DIR / f"{cat}.tsv"
        if tsv_path.exists():
            with open(tsv_path) as f:
                total_queries += sum(1 for _ in csv.DictReader(f, delimiter="\t"))
    total_queries_by_db[dbtype] = total_queries

# Generate statistics for each database
stats_by_db = {}
for dbtype in DB_TYPES:
    print(f"[STATS] Calculating statistics for {dbtype}...")
    stats = calculate_blast_statistics(all_hits_by_db[dbtype], dbtype, total_queries_by_db[dbtype])
    stats_by_db[dbtype] = stats
    write_blast_statistics(stats, OUT_DIR / dbtype, dbtype, category_stats_by_db[dbtype])

# Create comparison statistics
print("[STATS] Creating database comparison...")
create_comparison_statistics(stats_by_db['premature'], stats_by_db['mature'], OUT_DIR)

# -----------------------------------
# SUMMARY
# -----------------------------------
print("\n=== SUMMARY ===")
sumtxt = OUT_DIR / "summary.txt"
sumtsv = OUT_DIR / "summary.tsv"

with open(sumtxt, "w") as f, open(sumtsv, "w") as g:
    g.write("Class\tDB\tTotalHits\n")
    for cat, db, hits in summary:
        line = f"{cat}\t{db}\t{hits}"
        print(line)
        f.write(line + "\n")
        g.write(line + "\n")

print("\n=== BLAST_DANEEL_ULTRA v2 COMPLETE ===")
print("Enhanced statistics generated in:")
print("  - blast_results/premature/statistics/")
print("  - blast_results/mature/statistics/")
print("  - blast_results/comparison/")