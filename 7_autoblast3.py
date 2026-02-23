#!/usr/bin/env python3
"""
BLAST_DANEEL_ULTRA v5.4 - ENHANCED OPTIMIZATION!
====================================================================
★ E-value MIN and MAX are optimized separately!
★ Identity-Coverage combined score + correlation analysis!
★ MAXIMUM 5 hits per intron (highest bitscore)
★ Category-based realistic thresholds
★ TAKES 6_supervalidator_stat.py OUTPUT AS INPUT!
"""

import csv
import subprocess
import statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter, defaultdict
from Bio import SeqIO
import sys
import warnings
warnings.filterwarnings('ignore')

# ===========================================
# PLOTLY - optional
# ===========================================
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.io as pio
    pio.kaleido.scope.default_format = "png"
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

print("=" * 80)
print("BLAST_DANEEL_ULTRA v5.4 - ENHANCED OPTIMIZATION!")
print("=" * 80)
print("⚡ Maximum 5 hits per intron (highest bitscore)")
print("⚡ E-value MIN and MAX are optimized separately!")
print("⚡ Identity-Coverage combined score + correlation analysis!")
print("⚡ Takes 6_supervalidator_stat.py OUTPUT AS INPUT!")
print("=" * 80)

# ===========================================
# ASK USER - VALIDATOR OUTPUT
# ===========================================
print("\n📁 SPECIFY VALIDATOR OUTPUT FOLDER")
print("-" * 80)
print("When you run 6_supervalidator_stat.py, a 'validator_results/' folder is created.")
print("This folder should contain .fa and .tsv files under 'classes/'.")

default_validator_dir = "validator_results"
validator_input = input(f"\n➡️  Validator results folder [default: {default_validator_dir}]: ").strip()
VALIDATOR_DIR = Path(validator_input if validator_input else default_validator_dir)

CLASS_DIR = VALIDATOR_DIR / "classes"
if not CLASS_DIR.exists():
    print(f"\n❌ ERROR: {CLASS_DIR} folder not found!")
    print(f"   Did you run 6_supervalidator_stat.py first?")
    sys.exit(1)

# Check categories
for cat in ["very_high", "high", "moderate", "low"]:
    if not (CLASS_DIR / f"{cat}.fa").exists():
        print(f"\n❌ ERROR: {CLASS_DIR / f'{cat}.fa'} not found!")
        sys.exit(1)

print(f"\n✅ Validator output found: {VALIDATOR_DIR}")

# ===========================================
# ASK FOR DATABASE FASTA
# ===========================================
print(f"\n🎯 Database FASTA files?")
db_pre_input = input("➡️  premature.fasta path [default: premature.fasta]: ").strip()
db_mat_input = input("➡️  mature.fasta path [default: mature.fasta]: ").strip()

DB_PREMATURE = Path(db_pre_input if db_pre_input else "premature.fasta")
DB_MATURE = Path(db_mat_input if db_mat_input else "mature.fasta")

if not DB_PREMATURE.exists() or not DB_MATURE.exists():
    print(f"❌ ERROR: Database files not found!")
    sys.exit(1)

# ===========================================
# ASK FOR OUTPUT FOLDER
# ===========================================
default_out = "blast_results_ultimate"
out_dir_input = input(f"\n📂 Output folder [default: {default_out}]: ").strip()
OUT_DIR = Path(out_dir_input if out_dir_input else default_out)

# ===========================================
# ASK FOR RUNNING MODE
# ===========================================
print(f"\n⚙️  Running mode?")
print(f"   1) Fully automatic (optimization first, then BLAST) - RECOMMENDED")
print(f"   2) Run optimization only, skip BLAST")
print(f"   3) Run BLAST only (enter manual thresholds)")
mode = input("➡️  Choice (1/2/3): ").strip() or "1"

# ===========================================
# ASK FOR MANUAL THRESHOLDS
# ===========================================
manual_pident = None
manual_coverage = None
manual_evalue_min = None
manual_evalue_max = None

if mode == "3":
    print(f"\n🎛️  Enter manual thresholds:")
    manual_pident = float(input("➡️  Percent Identity (≥): ").strip())
    manual_coverage = float(input("➡️  Coverage (≥): ").strip())
    manual_evalue_min = float(input("➡️  E-value MIN (≥, e.g., 1e-10): ").strip())
    manual_evalue_max = float(input("➡️  E-value MAX (≤, e.g., 1e-2): ").strip())

print("\n" + "=" * 80)
print("✅ INPUTS RECEIVED!")
print(f"📁 Validator output: {VALIDATOR_DIR}")
print(f"📁 Class files: {CLASS_DIR}")
print(f"🎯 Database: {DB_PREMATURE.name}, {DB_MATURE.name}")
print(f"📂 Output: {OUT_DIR}")
print(f"⚙️  Mode: {mode}")
print("=" * 80)

# ===========================================
# GLOBAL PARAMETERS
# ===========================================
CATS = ["very_high", "high", "moderate", "low"]
DB_TYPES = ["premature", "mature"]

# Directories
OPT_DIR = OUT_DIR / "optimization"
HTML_DIR = OUT_DIR / "interactive_plots"
DB_DIR = OUT_DIR / "db"
PLOT_DIR = OUT_DIR / "plots"
VALIDATION_DIR = OUT_DIR / "validation"

for d in [OUT_DIR, OPT_DIR, HTML_DIR, DB_DIR, PLOT_DIR, VALIDATION_DIR,
          OUT_DIR / "premature", OUT_DIR / "mature",
          OUT_DIR / "premature/statistics", OUT_DIR / "mature/statistics",
          OUT_DIR / "premature/fasta", OUT_DIR / "mature/fasta"]:
    d.mkdir(parents=True, exist_ok=True)

HEADER_BLAST = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]

# ===========================================
# SIMPLE STATISTICS FUNCTIONS
# ===========================================
def percentile(data, p):
    """Simple percentile calculation"""
    if not data:
        return 0
    sorted_data = sorted(data)
    k = (len(sorted_data) - 1) * (p / 100)
    f = int(k)
    c = int(np.ceil(k))
    if f == c:
        return sorted_data[f]
    d0 = sorted_data[f] * (c - k)
    d1 = sorted_data[c] * (k - f)
    return d0 + d1

def elbow_point(x, y):
    """Find elbow point - largest second derivative"""
    if len(y) < 3:
        return len(y) // 2
    y = np.array(y)
    second_derivative = np.gradient(np.gradient(y))
    return np.argmax(np.abs(second_derivative))

# ===========================================
# NEW: CORRELATION ANALYSIS
# ===========================================
def correlation_analysis(pidents, coverages):
    """Calculate Pearson correlation between Identity and Coverage"""
    if len(pidents) < 3 or len(coverages) < 3:
        return 0, "Insufficient data"
    
    try:
        corr = np.corrcoef(pidents, coverages)[0, 1]
        
        if np.isnan(corr):
            return 0, "Correlation could not be calculated"
        
        if abs(corr) > 0.7:
            comment = f"⚠️ HIGH CORRELATION (|r|={corr:.3f}) - A single criterion might be sufficient"
        elif abs(corr) > 0.5:
            comment = f"📊 MODERATE CORRELATION (|r|={corr:.3f}) - Both criteria can be used together"
        else:
            comment = f"✅ LOW CORRELATION (|r|={corr:.3f}) - Both criteria provide independent information"
        
        return corr, comment
    except:
        return 0, "Correlation calculation error"

# ===========================================
# NEW: COMBINED SCORE
# ===========================================
def combined_score(pident, coverage, w1=0.7, w2=0.3):
    """
    Weighted score combining Identity and Coverage
    Default: Identity 70%, Coverage 30% weight
    """
    # Normalize (to 0-1 range)
    pident_norm = pident / 100.0
    coverage_norm = coverage / 100.0
    
    # Weighted combined score
    return w1 * pident_norm + w2 * coverage_norm

def optimize_combined_threshold(features, weights=(0.7, 0.3)):
    """
    Find optimal threshold for combined score
    """
    if len(features) < 10:
        return 0.7  # Default
    
    combined_scores = []
    for f in features:
        score = combined_score(f['pident'], f['coverage'], weights[0], weights[1])
        combined_scores.append(score)
    
    # Use 15th percentile (tolerant)
    combined_th = percentile(combined_scores, 15)
    
    return combined_th

# ===========================================
# UPDATED: E-VALUE RANGE OPTIMIZATION
# ===========================================
def optimize_evalue_range(evalues):
    """
    Find optimal MIN and MAX values for E-value.
    - MIN: 10th percentile (lowest 10%)
    - MAX: 90th percentile (highest 10%)
    
    """
    if not evalues:
        return 1e-10, 1e-2
    
    # MIN value - 10th percentile (lowest 10%)
    evalue_min = percentile(evalues, 10)
    
    # MAX value - 90th percentile (highest 10%)
    evalue_max = percentile(evalues, 90)
    
    # Safety limits
    if evalue_min < 1e-50:  # Too small
        evalue_min = 1e-50
    if evalue_min > 1e-5:    # Too large
        evalue_min = 1e-5
        
    if evalue_max > 0.1:      # Too large
        evalue_max = 0.1
    if evalue_max < 1e-6:     # Too small
        evalue_max = 1e-6
    
    return evalue_min, evalue_max

def optimize_evalue_range_category(evalues):
    """
    Category-based E-value range optimization.
    - MIN: 5th percentile (more tolerant)
    - MAX: 85th percentile
    """
    if not evalues:
        return 1e-12, 1e-3
    
    evalue_min = percentile(evalues, 5)   # 5th percentile (more tolerant)
    evalue_max = percentile(evalues, 85)  # 85th percentile
    
    if evalue_min < 1e-50:
        evalue_min = 1e-50
    if evalue_min > 1e-5:
        evalue_min = 1e-5
        
    if evalue_max > 0.1:
        evalue_max = 0.1
    if evalue_max < 1e-7:
        evalue_max = 1e-7
    
    return evalue_min, evalue_max

# ===========================================
# MULTIPLE MATCH FILTER
# ===========================================
def filter_multiple_matches(hits, max_matches=5):
    """
    Select at most 5 hits per intron (query).
    Selection criteria: BITSCORE (highest)
    """
    if not hits:
        return []
    
    query_groups = defaultdict(list)
    for hit in hits:
        query_id = hit[0]
        query_groups[query_id].append(hit)
    
    filtered_hits = []
    for query_id, group_hits in query_groups.items():
        sorted_hits = sorted(group_hits, key=lambda x: float(x[11]), reverse=True)
        best_hits = sorted_hits[:max_matches]
        filtered_hits.extend(best_hits)
    
    return filtered_hits

# ===========================================
# 1. EXPLORATORY BLAST - FINE TUNED!
# ===========================================
class ExploratoryBLAST:
    @staticmethod
    def run():
        print("\n[1/8] 🔍 PERFORMING BROAD SCAN...")
        print("-" * 80)
        print(f"   • Percent Identity: ≥ 70%")
        print(f"   • E-value: ≤ 0.01")
        print(f"   • Coverage filter: NONE")
        print(f"   • Query: {CLASS_DIR}/[cat].fa")
        print(f"   • Database: {DB_PREMATURE.name}, {DB_MATURE.name}")
        print("-" * 80)
        
        # Set up databases
        print("\n   📦 Setting up databases...")
        ExploratoryBLAST._build_db(str(DB_PREMATURE), DB_DIR / "db_premature")
        ExploratoryBLAST._build_db(str(DB_MATURE), DB_DIR / "db_mature")
        
        all_hits = []
        
        for cat in CATS:
            fasta = CLASS_DIR / f"{cat}.fa"
            tsv = CLASS_DIR / f"{cat}.tsv"
            
            if not fasta.exists() or not tsv.exists():
                continue
            
            print(f"\n   📁 Processing {cat}...")
            seqs = ExploratoryBLAST._parse_fasta(fasta)
            
            with open(tsv) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            
            for dbtype in DB_TYPES:
                hits = ExploratoryBLAST._run_blast(cat, dbtype, seqs, rows)
                all_hits.extend(hits)
                print(f"     • {dbtype}: {len(hits)} hits")
        
        print(f"\n   ✅ Total raw hits: {len(all_hits)}")
        return all_hits
    
    @staticmethod
    def _build_db(fasta, dbname):
        try:
            subprocess.run([
                "makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", str(dbname)
            ], check=True, capture_output=True, text=True)
            print(f"     ✅ Database created: {dbname}")
            return True
        except Exception as e:
            print(f"     ❌ Database error: {e}")
            return False
    
    @staticmethod
    def _parse_fasta(fapath):
        seqmap = {}
        try:
            for rec in SeqIO.parse(fapath, "fasta"):
                header = rec.id
                coord = header.split("|")[0] if "|" in header else header
                seqmap.setdefault(coord, []).append(str(rec.seq))
        except:
            pass
        return seqmap
    
    @staticmethod
    def _run_blast(cat, dbtype, seqs, rows):
        tempfa = Path(f"temp_explore_{cat}_{dbtype}.fa")
        try:
            with open(tempfa, "w") as tf:
                for r in rows:
                    cid = r["CoordID"]
                    if cid in seqs:
                        for s in seqs[cid][:1]:
                            tf.write(f">{cid}\n{s}\n")
            
            cmd = [
                "blastn", "-task", "blastn-short",
                "-query", str(tempfa),
                "-db", str(DB_DIR / f"db_{dbtype}"),
                "-outfmt", "6 " + " ".join(HEADER_BLAST),
                "-perc_identity", "70",
                "-evalue", "0.01"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            lines = result.stdout.strip().split("\n")
            
            hits = []
            if lines and lines[0]:
                for line in lines:
                    parts = line.split("\t")
                    if len(parts) >= 12:
                        hits.append({
                            'line': line,
                            'cat': cat,
                            'db': dbtype
                        })
            return hits
        except Exception as e:
            print(f"     ⚠️ BLAST error: {e}")
            return []
        finally:
            if tempfa.exists():
                tempfa.unlink()

# ===========================================
# 2. FEATURE EXTRACTOR - UPDATED
# ===========================================
class FeatureExtractor:
    @staticmethod
    def extract(exploratory_hits):
        print("\n[2/8] 📊 ANALYZING DATA...")
        print("-" * 80)
        
        all_features = []
        category_features = {cat: [] for cat in CATS}
        
        for hit_data in exploratory_hits:
            hit = hit_data['line'].split("\t")
            try:
                pident = float(hit[2])
                length = int(hit[3])
                qstart = int(hit[6])
                qend = int(hit[7])
                query_length = abs(qend - qstart) + 1
                coverage = (length / query_length) * 100 if query_length > 0 else 0
                evalue = float(hit[10])
                bitscore = float(hit[11])
                
                if 0 < coverage <= 100 and 0 < pident <= 100 and evalue > 0:
                    # Calculate combined score
                    combined = combined_score(pident, coverage)
                    
                    feature = {
                        'pident': pident,
                        'coverage': coverage,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'combined_score': combined,
                        'cat': hit_data['cat'],
                        'db': hit_data['db']
                    }
                    all_features.append(feature)
                    if hit_data['cat'] in category_features:
                        category_features[hit_data['cat']].append(feature)
            except Exception as e:
                continue
        
        print(f"   • Valid data points: {len(all_features)}")
        if all_features:
            print(f"   • Average pident: {np.mean([f['pident'] for f in all_features]):.1f}%")
            print(f"   • Average coverage: {np.mean([f['coverage'] for f in all_features]):.1f}%")
            print(f"   • Average combined score: {np.mean([f['combined_score'] for f in all_features]):.3f}")
            print(f"   • Average evalue: {np.mean([f['evalue'] for f in all_features]):.2e}")
            
            # Correlation analysis
            pidents = [f['pident'] for f in all_features]
            coverages = [f['coverage'] for f in all_features]
            corr, comment = correlation_analysis(pidents, coverages)
            print(f"   • Identity-Coverage correlation: {corr:.3f}")
            print(f"   • {comment}")
        
        return all_features, category_features

# ===========================================
# 3. THRESHOLD OPTIMIZER - UPDATED
# ===========================================
class ThresholdOptimizer:
    @staticmethod
    def optimize(features):
        """Find threshold using Percentile + Elbow - COMBINED SCORE SUPPORT"""
        
        if len(features) < 10:
            return None
        
        pidents = [f['pident'] for f in features]
        coverages = [f['coverage'] for f in features]
        evalues = [f['evalue'] for f in features]
        combined_scores = [f['combined_score'] for f in features]
        
        # Correlation analysis
        corr, _ = correlation_analysis(pidents, coverages)
        
        # Method 1: Percentile
        pident_pctl = percentile(pidents, 15)
        coverage_pctl = percentile(coverages, 15)
        combined_pctl = percentile(combined_scores, 15)
        
        # Method 2: Elbow (for pident)
        pident_range = range(60, 100, 2)
        hit_counts = []
        for th in pident_range:
            count = sum(1 for f in features if f['pident'] >= th)
            hit_counts.append(count)
        
        elbow_idx = elbow_point(list(pident_range), hit_counts)
        pident_elbow = pident_range[elbow_idx]
        
        # Ensemble: average of both
        pident_final = (pident_pctl + pident_elbow) / 2
        
        # Optimize E-value range - DATA-DRIVEN!
        evalue_min, evalue_max = optimize_evalue_range(evalues)
        
        # If high correlation, recommend combined score
        if abs(corr) > 0.7:
            method_note = f"ENSEMBLE + HIGH CORRELATION (r={corr:.2f}) - COMBINED SCORE RECOMMENDED: {combined_pctl:.3f}"
        else:
            method_note = f"ENSEMBLE + LOW CORRELATION (r={corr:.2f})"
        
        return {
            'pident': pident_final,
            'coverage': coverage_pctl,
            'combined_score': combined_pctl,
            'evalue_min': evalue_min,
            'evalue_max': evalue_max,
            'correlation': corr,
            'method': method_note
        }
    
    @staticmethod
    def optimize_category(features):
        """Category-based threshold - COMBINED SCORE SUPPORT"""
        if len(features) < 5:
            return None
        
        pidents = [f['pident'] for f in features]
        coverages = [f['coverage'] for f in features]
        evalues = [f['evalue'] for f in features]
        combined_scores = [f['combined_score'] for f in features]
        
        # Correlation analysis
        corr, _ = correlation_analysis(pidents, coverages)
        
        # More tolerant for category
        pident_th = percentile(pidents, 10)
        coverage_th = percentile(coverages, 10)
        combined_th = percentile(combined_scores, 10)
        
        # Optimize E-value range - category-based, DATA-DRIVEN!
        evalue_min, evalue_max = optimize_evalue_range_category(evalues)
        
        return {
            'pident': pident_th,
            'coverage': coverage_th,
            'combined_score': combined_th,
            'evalue_min': evalue_min,
            'evalue_max': evalue_max,
            'correlation': corr
        }

# ===========================================
# 4. CROSS VALIDATOR
# ===========================================
class CrossValidator:
    @staticmethod
    def validate(features):
        print("\n[3/8] 🔬 PERFORMING CROSS-VALIDATION...")
        print("-" * 80)
        
        if len(features) < 50:
            print("   ⚠️ Insufficient data, skipping cross-validation")
            return None
        
        n_folds = 5
        fold_size = len(features) // n_folds
        pident_folds = []
        
        for i in range(n_folds):
            start = i * fold_size
            end = start + fold_size if i < n_folds - 1 else len(features)
            fold_features = features[start:end]
            pidents = [f['pident'] for f in fold_features]
            pident_folds.append(percentile(pidents, 15))
        
        return {
            'pident_mean': np.mean(pident_folds),
            'pident_std': np.std(pident_folds),
            'stability_score': 1 / (1 + np.std(pident_folds))
        }

# ===========================================
# 5. PLOT GENERATOR - UPDATED
# ===========================================
class PlotGenerator:
    @staticmethod
    def generate(features, thresholds, output_dir):
        print("\n[4/8] 📈 PREPARING PLOTS...")
        
        plt.figure(figsize=(18, 6))
        
        # 1. Identity vs Coverage (with correlation)
        plt.subplot(1, 3, 1)
        categories = [f['cat'] for f in features]
        pidents = [f['pident'] for f in features]
        coverages = [f['coverage'] for f in features]
        
        colors = {'very_high': 'darkgreen', 'high': 'steelblue', 
                  'moderate': 'orange', 'low': 'red'}
        cat_colors = [colors.get(cat, 'gray') for cat in categories]
        
        plt.scatter(pidents, coverages, c=cat_colors, alpha=0.6, s=30)
        plt.axvline(x=thresholds['pident'], color='red', linestyle='--', 
                   label=f"pident ≥ {thresholds['pident']:.1f}%")
        plt.axhline(y=thresholds['coverage'], color='red', linestyle='--',
                   label=f"coverage ≥ {thresholds['coverage']:.1f}%")
        
        # Add correlation info
        corr = thresholds.get('correlation', 0)
        plt.text(0.05, 0.95, f"Correlation: r={corr:.3f}", 
                transform=plt.gca().transAxes, fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plt.xlabel('Percent Identity (%)')
        plt.ylabel('Coverage (%)')
        plt.title('Identity vs Coverage')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 2. Combined Score Distribution
        plt.subplot(1, 3, 2)
        combined_scores = [f.get('combined_score', 0) for f in features]
        
        plt.hist(combined_scores, bins=50, alpha=0.7, color='purple', edgecolor='black')
        plt.axvline(x=thresholds.get('combined_score', 0.7), color='red', linestyle='--',
                   label=f"Combined ≥ {thresholds.get('combined_score', 0.7):.3f}")
        plt.xlabel('Combined Score (0.7*Identity + 0.3*Coverage)')
        plt.ylabel('Frequency')
        plt.title('Combined Score Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 3. E-value Distribution
        plt.subplot(1, 3, 3)
        evalues = [f['evalue'] for f in features]
        log_evalues = [-np.log10(e + 1e-300) for e in evalues]
        
        plt.hist(log_evalues, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        plt.axvline(x=-np.log10(thresholds['evalue_min']), color='green', linestyle='--',
                   label=f"MIN: {thresholds['evalue_min']:.1e}")
        plt.axvline(x=-np.log10(thresholds['evalue_max']), color='red', linestyle='--',
                   label=f"MAX: {thresholds['evalue_max']:.1e}")
        plt.xlabel('-log10(E-value)')
        plt.ylabel('Frequency')
        plt.title('E-value Distribution (Data-Driven Range)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'threshold_optimization.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"   ✅ Plot saved: {output_dir}/threshold_optimization.png")

# ===========================================
# 6. FINAL BLAST - COMBINED SCORE SUPPORT
# ===========================================
class FinalBLAST:
    @staticmethod
    def run(thresholds, category_thresholds=None):
        print("\n[5/8] 🚀 STARTING FINAL BLAST...")
        print(f"   • Global pident ≥ {thresholds['pident']:.1f}%")
        print(f"   • Global coverage ≥ {thresholds['coverage']:.1f}%")
        print(f"   • Global combined score ≥ {thresholds.get('combined_score', 0.7):.3f}")
        print(f"   • Global evalue: {thresholds['evalue_min']:.1e} - {thresholds['evalue_max']:.1e} (Data-Driven!)")
        print(f"   • Maximum hits/intron: 5 (highest bitscore)")
        print("-" * 80)
        
        summary = []
        total_hits_all = 0
        total_introns_all = 0
        
        for cat in CATS:
            fasta = CLASS_DIR / f"{cat}.fa"
            tsv = CLASS_DIR / f"{cat}.tsv"
            
            if not fasta.exists() or not tsv.exists():
                continue
            
            print(f"\n   📁 Processing {cat}...")
            seqs = ExploratoryBLAST._parse_fasta(fasta)
            
            with open(tsv) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            
            total_introns = len(rows)
            total_introns_all += total_introns
            
            # Category-specific threshold
            if category_thresholds and cat in category_thresholds:
                cat_th = category_thresholds[cat]
                pident_th = cat_th['pident']
                coverage_th = cat_th['coverage']
                evalue_min_th = cat_th['evalue_min']
                evalue_max_th = cat_th['evalue_max']
                print(f"     • Category-specific: pident≥{pident_th:.1f}%, cov≥{coverage_th:.1f}%, "
                      f"evalue: {evalue_min_th:.1e}-{evalue_max_th:.1e}")
            else:
                pident_th = thresholds['pident']
                coverage_th = thresholds['coverage']
                evalue_min_th = thresholds['evalue_min']
                evalue_max_th = thresholds['evalue_max']
            
            for dbtype in DB_TYPES:
                hits = FinalBLAST._run_blast_category(
                    cat, dbtype, seqs, rows,
                    pident_th, coverage_th, evalue_min_th, evalue_max_th
                )
                
                hits_after_filter = len(hits)
                total_hits_all += hits_after_filter
                
                avg_hits_per_intron = hits_after_filter / total_introns if total_introns > 0 else 0
                summary.append((cat, dbtype, hits_after_filter, total_introns))
                print(f"     • {dbtype}: {hits_after_filter} hits "
                      f"({avg_hits_per_intron:.2f} avg hits/intron)")
        
        print(f"\n   📊 OVERALL STATISTICS:")
        print(f"     • Total introns: {total_introns_all}")
        print(f"     • Total hits: {total_hits_all}")
        print(f"     • Average hits/intron: {total_hits_all/total_introns_all:.2f}")
        print(f"     • Maximum hits/intron: 5 (rule applied)")
        
        return summary
    
    @staticmethod
    def _run_blast_category(cat, dbtype, seqs, rows, pident_th, coverage_th, evalue_min_th, evalue_max_th):
        tempfa = Path(f"temp_final_{cat}_{dbtype}.fa")
        try:
            with open(tempfa, "w") as tf:
                for r in rows:
                    cid = r["CoordID"]
                    if cid in seqs:
                        for s in seqs[cid][:1]:
                            tf.write(f">{cid}\n{s}\n")
            
            cmd = [
                "blastn", "-task", "blastn-short",
                "-query", str(tempfa),
                "-db", str(DB_DIR / f"db_{dbtype}"),
                "-outfmt", "6 " + " ".join(HEADER_BLAST),
                "-perc_identity", str(pident_th),
                "-evalue", str(evalue_max_th)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            lines = result.stdout.strip().split("\n")
            raw_hits = lines if lines and lines[0] else []
            
            # First filter by evalue range (between MIN and MAX)
            evalue_filtered = []
            for line in raw_hits:
                parts = line.split("\t")
                if len(parts) < 12:
                    continue
                try:
                    evalue = float(parts[10])
                    if evalue_min_th <= evalue <= evalue_max_th:
                        evalue_filtered.append(parts)
                except:
                    continue
            
            # Then filter by coverage
            coverage_filtered = []
            for parts in evalue_filtered:
                length = int(parts[3])
                qstart = int(parts[6])
                qend = int(parts[7])
                query_length = abs(qend - qstart) + 1
                coverage = (length / query_length) * 100 if query_length > 0 else 0
                
                if coverage >= coverage_th:
                    coverage_filtered.append(parts)
            
            # Multiple match filter
            filtered_hits = filter_multiple_matches(coverage_filtered, max_matches=5)
            
            # Write TSV
            out_tsv = OUT_DIR / dbtype / f"{cat}.tsv"
            out_tsv.parent.mkdir(exist_ok=True)
            
            with open(out_tsv, "w", newline="") as f:
                if rows:
                    fn = list(rows[0].keys()) + [
                        "BLAST_HIT", "pident", "evalue", "bitscore",
                        "qstart", "qend", "sstart", "send", "DB", "coverage", "combined_score"
                    ]
                    w = csv.DictWriter(f, fieldnames=fn, delimiter="\t")
                    w.writeheader()
                    
                    hit_map = defaultdict(list)
                    for h in filtered_hits:
                        hit_map[h[0]].append(h)
                    
                    for r in rows:
                        cid = r["CoordID"]
                        if cid in hit_map:
                            for h in hit_map[cid]:
                                new = r.copy()
                                _, sid, pid, ln, mis, gap, qs, qe, ss, se, ev, bit = h
                                qlen = abs(int(qe) - int(qs)) + 1
                                cov = (int(ln) / qlen) * 100
                                combined = combined_score(float(pid), cov)
                                new.update({
                                    "BLAST_HIT": sid,
                                    "pident": pid,
                                    "evalue": ev,
                                    "bitscore": bit,
                                    "qstart": qs,
                                    "qend": qe,
                                    "sstart": ss,
                                    "send": se,
                                    "DB": dbtype,
                                    "coverage": f"{cov:.1f}",
                                    "combined_score": f"{combined:.3f}"
                                })
                                w.writerow(new)
                        else:
                            new = r.copy()
                            new.update({
                                "BLAST_HIT": "", "pident": "", "evalue": "", "bitscore": "",
                                "qstart": "", "qend": "", "sstart": "", "send": "", 
                                "DB": dbtype, "coverage": "", "combined_score": ""
                            })
                            w.writerow(new)
            
            # Write FASTA
            out_fa = OUT_DIR / dbtype / "fasta" / f"{cat}.fa"
            out_fa.parent.mkdir(parents=True, exist