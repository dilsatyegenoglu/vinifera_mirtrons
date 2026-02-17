#!/usr/bin/env python3
"""
BLAST_DANEEL_ULTRA v5.4 - ENHANCED OPTIMIZATION!
====================================================================
★ E-value MIN ve MAX ayrı ayrı optimize edilir (veriye dayalı)!
★ Identity-Coverage kombine skor + korelasyon analizi!
★ Her intron için EN FAZLA 5 hit (en yüksek bitscore)
★ Kategori bazlı gerçekçi threshold'lar
★ 6_supervalidator_stat.py OUTPUT'UNU INPUT ALIR!
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
# PLOTLY - opsiyonel
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
print("⚡ Her intron için EN FAZLA 5 hit (en yüksek bitscore)")
print("⚡ E-value MIN ve MAX ayrı ayrı optimize edilir (veriye dayalı)!")
print("⚡ Identity-Coverage kombine skor + korelasyon analizi!")
print("⚡ 6_supervalidator_stat.py OUTPUT'UNU INPUT ALIR!")
print("=" * 80)

# ===========================================
# KULLANICIYA SOR - VALIDATOR OUTPUT
# ===========================================
print("\n📁 VALIDATOR OUTPUT KLASÖRÜNÜ BELİRLEYİN")
print("-" * 80)
print("6_supervalidator_stat.py çalıştırdığınızda 'validator_results/' klasörü oluşur.")
print("Bu klasörde 'classes/' altında .fa ve .tsv dosyaları olmalı.")

default_validator_dir = "validator_results"
validator_input = input(f"\n➡️  Validator results klasörü [varsayılan: {default_validator_dir}]: ").strip()
VALIDATOR_DIR = Path(validator_input if validator_input else default_validator_dir)

CLASS_DIR = VALIDATOR_DIR / "classes"
if not CLASS_DIR.exists():
    print(f"\n❌ HATA: {CLASS_DIR} klasörü bulunamadı!")
    print(f"   Önce 6_supervalidator_stat.py çalıştırdınız mı?")
    sys.exit(1)

# Kategorileri kontrol et
for cat in ["very_high", "high", "moderate", "low"]:
    if not (CLASS_DIR / f"{cat}.fa").exists():
        print(f"\n❌ HATA: {CLASS_DIR / f'{cat}.fa'} bulunamadı!")
        sys.exit(1)

print(f"\n✅ Validator output bulundu: {VALIDATOR_DIR}")

# ===========================================
# DATABASE FASTA SOR
# ===========================================
print(f"\n🎯 Database FASTA dosyaları?")
db_pre_input = input("➡️  premature.fasta yolu [varsayılan: premature.fasta]: ").strip()
db_mat_input = input("➡️  mature.fasta yolu [varsayılan: mature.fasta]: ").strip()

DB_PREMATURE = Path(db_pre_input if db_pre_input else "premature.fasta")
DB_MATURE = Path(db_mat_input if db_mat_input else "mature.fasta")

if not DB_PREMATURE.exists() or not DB_MATURE.exists():
    print(f"❌ HATA: Database dosyaları bulunamadı!")
    sys.exit(1)

# ===========================================
# OUTPUT KLASÖRÜ SOR
# ===========================================
default_out = "blast_results_ultimate"
out_dir_input = input(f"\n📂 Output klasörü [varsayılan: {default_out}]: ").strip()
OUT_DIR = Path(out_dir_input if out_dir_input else default_out)

# ===========================================
# ÇALIŞMA MODU SOR
# ===========================================
print(f"\n⚙️  Çalışma modu?")
print(f"   1) Tam otomatik (önce optimizasyon, sonra BLAST) - ÖNERİLEN")
print(f"   2) Sadece optimizasyon yap, BLAST yapma")
print(f"   3) Sadece BLAST yap (manuel eşik gir)")
mode = input("➡️  Seçim (1/2/3): ").strip() or "1"

# ===========================================
# MANUEL MODDA EŞİK SOR
# ===========================================
manual_pident = None
manual_coverage = None
manual_evalue_min = None
manual_evalue_max = None

if mode == "3":
    print(f"\n🎛️  Manuel eşikleri girin:")
    manual_pident = float(input("➡️  Percent Identity (≥): ").strip())
    manual_coverage = float(input("➡️  Coverage (≥): ").strip())
    manual_evalue_min = float(input("➡️  E-value MIN (≥, örn: 1e-10): ").strip())
    manual_evalue_max = float(input("➡️  E-value MAX (≤, örn: 1e-2): ").strip())

print("\n" + "=" * 80)
print("✅ INPUTLAR ALINDI!")
print(f"📁 Validator output: {VALIDATOR_DIR}")
print(f"📁 Class files: {CLASS_DIR}")
print(f"🎯 Database: {DB_PREMATURE.name}, {DB_MATURE.name}")
print(f"📂 Output: {OUT_DIR}")
print(f"⚙️  Mod: {mode}")
print("=" * 80)

# ===========================================
# GLOBAL PARAMETRELER
# ===========================================
CATS = ["very_high", "high", "moderate", "low"]
DB_TYPES = ["premature", "mature"]

# Dizinler
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
# BASİT İSTATİSTİK FONKSİYONLARI
# ===========================================
def percentile(data, p):
    """Basit percentile hesaplama"""
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
    """Elbow point bul - en büyük ikinci türev"""
    if len(y) < 3:
        return len(y) // 2
    y = np.array(y)
    second_derivative = np.gradient(np.gradient(y))
    return np.argmax(np.abs(second_derivative))

# ===========================================
# YENİ: KORELASYON ANALİZİ
# ===========================================
def correlation_analysis(pidents, coverages):
    """Identity ile Coverage arasındaki Pearson korelasyonunu hesapla"""
    if len(pidents) < 3 or len(coverages) < 3:
        return 0, "Yetersiz veri"
    
    try:
        corr = np.corrcoef(pidents, coverages)[0, 1]
        
        if np.isnan(corr):
            return 0, "Korelasyon hesaplanamadı"
        
        if abs(corr) > 0.7:
            yorum = f"⚠️ YÜKSEK KORELASYON (|r|={corr:.3f}) - Tek kriter yeterli olabilir"
        elif abs(corr) > 0.5:
            yorum = f"📊 ORTA KORELASYON (|r|={corr:.3f}) - İki kriter birlikte kullanılabilir"
        else:
            yorum = f"✅ DÜŞÜK KORELASYON (|r|={corr:.3f}) - İki kriter bağımsız bilgi veriyor"
        
        return corr, yorum
    except:
        return 0, "Korelasyon hesaplama hatası"

# ===========================================
# YENİ: KOMBİNE SKOR
# ===========================================
def combined_score(pident, coverage, w1=0.7, w2=0.3):
    """
    Identity ve Coverage'ı birleştiren weighted skor
    Varsayılan: Identity %70, Coverage %30 ağırlık
    """
    # Normalize et (0-1 aralığına)
    pident_norm = pident / 100.0
    coverage_norm = coverage / 100.0
    
    # Weighted combined score
    return w1 * pident_norm + w2 * coverage_norm

def optimize_combined_threshold(features, weights=(0.7, 0.3)):
    """
    Combined skor için optimal threshold bul
    """
    if len(features) < 10:
        return 0.7  # Varsayılan
    
    combined_scores = []
    for f in features:
        score = combined_score(f['pident'], f['coverage'], weights[0], weights[1])
        combined_scores.append(score)
    
    # 15. percentile kullan (toleranslı)
    combined_th = percentile(combined_scores, 15)
    
    return combined_th

# ===========================================
# GÜNCELLENMİŞ: E-VALUE RANGE OPTIMIZATION
# ===========================================
def optimize_evalue_range(evalues):
    """
    E-value için optimal MIN ve MAX değerlerini bul.
    - MIN: 10. percentile (en düşük %10)
    - MAX: 90. percentile (en yüksek %10)
    - Veriye dayalı optimizasyon!
    """
    if not evalues:
        return 1e-10, 1e-2
    
    # MIN değer - 10. percentile (en düşük %10)
    evalue_min = percentile(evalues, 10)
    
    # MAX değer - 90. percentile (en yüksek %10)
    evalue_max = percentile(evalues, 90)
    
    # Güvenlik sınırlamaları
    if evalue_min < 1e-50:  # Çok küçük çıkarsa
        evalue_min = 1e-50
    if evalue_min > 1e-5:    # Çok büyük çıkarsa
        evalue_min = 1e-5
        
    if evalue_max > 0.1:      # Çok büyük çıkarsa
        evalue_max = 0.1
    if evalue_max < 1e-6:     # Çok küçük çıkarsa
        evalue_max = 1e-6
    
    return evalue_min, evalue_max

def optimize_evalue_range_category(evalues):
    """
    Kategori bazlı E-value range optimizasyonu.
    - MIN: 5. percentile (daha toleranslı)
    - MAX: 85. percentile
    """
    if not evalues:
        return 1e-12, 1e-3
    
    evalue_min = percentile(evalues, 5)   # 5. percentile (daha toleranslı)
    evalue_max = percentile(evalues, 85)  # 85. percentile
    
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
    Her intron (query) için en fazla 5 hit seç.
    Seçim kriteri: BITSCORE (en yüksek)
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
        print("\n[1/8] 🔍 GENİŞ TARAMA YAPILIYOR...")
        print("-" * 80)
        print(f"   • Percent Identity: ≥ 70%")
        print(f"   • E-value: ≤ 0.01")
        print(f"   • Coverage filtresi: YOK")
        print(f"   • Query: {CLASS_DIR}/[cat].fa")
        print(f"   • Database: {DB_PREMATURE.name}, {DB_MATURE.name}")
        print("-" * 80)
        
        # Databaseleri kur
        print("\n   📦 Database'ler kuruluyor...")
        ExploratoryBLAST._build_db(str(DB_PREMATURE), DB_DIR / "db_premature")
        ExploratoryBLAST._build_db(str(DB_MATURE), DB_DIR / "db_mature")
        
        all_hits = []
        
        for cat in CATS:
            fasta = CLASS_DIR / f"{cat}.fa"
            tsv = CLASS_DIR / f"{cat}.tsv"
            
            if not fasta.exists() or not tsv.exists():
                continue
            
            print(f"\n   📁 {cat} işleniyor...")
            seqs = ExploratoryBLAST._parse_fasta(fasta)
            
            with open(tsv) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            
            for dbtype in DB_TYPES:
                hits = ExploratoryBLAST._run_blast(cat, dbtype, seqs, rows)
                all_hits.extend(hits)
                print(f"     • {dbtype}: {len(hits)} hit")
        
        print(f"\n   ✅ Toplam ham hit: {len(all_hits)}")
        return all_hits
    
    @staticmethod
    def _build_db(fasta, dbname):
        try:
            subprocess.run([
                "makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", str(dbname)
            ], check=True, capture_output=True, text=True)
            print(f"     ✅ Database oluşturuldu: {dbname}")
            return True
        except Exception as e:
            print(f"     ❌ Database hatası: {e}")
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
            print(f"     ⚠️ BLAST hatası: {e}")
            return []
        finally:
            if tempfa.exists():
                tempfa.unlink()

# ===========================================
# 2. FEATURE EXTRACTOR - GÜNCELLENDİ
# ===========================================
class FeatureExtractor:
    @staticmethod
    def extract(exploratory_hits):
        print("\n[2/8] 📊 VERİ ANALİZ EDİLİYOR...")
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
                    # Combined score hesapla
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
        
        print(f"   • Geçerli veri noktası: {len(all_features)}")
        if all_features:
            print(f"   • Ortalama pident: {np.mean([f['pident'] for f in all_features]):.1f}%")
            print(f"   • Ortalama coverage: {np.mean([f['coverage'] for f in all_features]):.1f}%")
            print(f"   • Ortalama combined score: {np.mean([f['combined_score'] for f in all_features]):.3f}")
            print(f"   • Ortalama evalue: {np.mean([f['evalue'] for f in all_features]):.2e}")
            
            # Korelasyon analizi
            pidents = [f['pident'] for f in all_features]
            coverages = [f['coverage'] for f in all_features]
            corr, yorum = correlation_analysis(pidents, coverages)
            print(f"   • Identity-Coverage korelasyonu: {corr:.3f}")
            print(f"   • {yorum}")
        
        return all_features, category_features

# ===========================================
# 3. THRESHOLD OPTIMIZER - GÜNCELLENDİ
# ===========================================
class ThresholdOptimizer:
    @staticmethod
    def optimize(features):
        """Percentile + Elbow ile threshold bul - KOMBİNE SKOR DESTEKLİ"""
        
        if len(features) < 10:
            return None
        
        pidents = [f['pident'] for f in features]
        coverages = [f['coverage'] for f in features]
        evalues = [f['evalue'] for f in features]
        combined_scores = [f['combined_score'] for f in features]
        
        # Korelasyon analizi
        corr, _ = correlation_analysis(pidents, coverages)
        
        # Method 1: Percentile
        pident_pctl = percentile(pidents, 15)
        coverage_pctl = percentile(coverages, 15)
        combined_pctl = percentile(combined_scores, 15)
        
        # Method 2: Elbow (pident için)
        pident_range = range(60, 100, 2)
        hit_counts = []
        for th in pident_range:
            count = sum(1 for f in features if f['pident'] >= th)
            hit_counts.append(count)
        
        elbow_idx = elbow_point(list(pident_range), hit_counts)
        pident_elbow = pident_range[elbow_idx]
        
        # Ensemble: ikisinin ortalaması
        pident_final = (pident_pctl + pident_elbow) / 2
        
        # E-value range optimize et - VERİYE DAYALI!
        evalue_min, evalue_max = optimize_evalue_range(evalues)
        
        # Yüksek korelasyon varsa combined score öner
        if abs(corr) > 0.7:
            method_note = f"ENSEMBLE + YÜKSEK KORELASYON (r={corr:.2f}) - COMBINED SKOR ÖNERİLİR: {combined_pctl:.3f}"
        else:
            method_note = f"ENSEMBLE + DÜŞÜK KORELASYON (r={corr:.2f})"
        
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
        """Kategori bazlı threshold - KOMBİNE SKOR DESTEKLİ"""
        if len(features) < 5:
            return None
        
        pidents = [f['pident'] for f in features]
        coverages = [f['coverage'] for f in features]
        evalues = [f['evalue'] for f in features]
        combined_scores = [f['combined_score'] for f in features]
        
        # Korelasyon analizi
        corr, _ = correlation_analysis(pidents, coverages)
        
        # Kategori için daha toleranslı
        pident_th = percentile(pidents, 10)
        coverage_th = percentile(coverages, 10)
        combined_th = percentile(combined_scores, 10)
        
        # E-value range optimize et - kategori bazlı, VERİYE DAYALI!
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
        print("\n[3/8] 🔬 CROSS-VALIDATION YAPILIYOR...")
        print("-" * 80)
        
        if len(features) < 50:
            print("   ⚠️ Yeterli veri yok, cross-validation atlanıyor")
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
# 5. PLOT GENERATOR - GÜNCELLENDİ
# ===========================================
class PlotGenerator:
    @staticmethod
    def generate(features, thresholds, output_dir):
        print("\n[4/8] 📈 GRAFİKLER HAZIRLANIYOR...")
        
        plt.figure(figsize=(18, 6))
        
        # 1. Identity vs Coverage (korelasyonlu)
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
        
        # Korelasyon bilgisini ekle
        corr = thresholds.get('correlation', 0)
        plt.text(0.05, 0.95, f"Korelasyon: r={corr:.3f}", 
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
        plt.title('E-value Distribution (Veriye Dayalı Range)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'threshold_optimization.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"   ✅ Grafik kaydedildi: {output_dir}/threshold_optimization.png")

# ===========================================
# 6. FINAL BLAST - KOMBİNE SKOR DESTEKLİ
# ===========================================
class FinalBLAST:
    @staticmethod
    def run(thresholds, category_thresholds=None):
        print("\n[5/8] 🚀 FİNAL BLAST BAŞLATILIYOR...")
        print(f"   • Global pident ≥ {thresholds['pident']:.1f}%")
        print(f"   • Global coverage ≥ {thresholds['coverage']:.1f}%")
        print(f"   • Global combined score ≥ {thresholds.get('combined_score', 0.7):.3f}")
        print(f"   • Global evalue: {thresholds['evalue_min']:.1e} - {thresholds['evalue_max']:.1e} (Veriye Dayalı!)")
        print(f"   • Maksimum hit/intron: 5 (en yüksek bitscore)")
        print("-" * 80)
        
        summary = []
        total_hits_all = 0
        total_introns_all = 0
        
        for cat in CATS:
            fasta = CLASS_DIR / f"{cat}.fa"
            tsv = CLASS_DIR / f"{cat}.tsv"
            
            if not fasta.exists() or not tsv.exists():
                continue
            
            print(f"\n   📁 {cat} işleniyor...")
            seqs = ExploratoryBLAST._parse_fasta(fasta)
            
            with open(tsv) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            
            total_introns = len(rows)
            total_introns_all += total_introns
            
            # Kategori bazlı threshold
            if category_thresholds and cat in category_thresholds:
                cat_th = category_thresholds[cat]
                pident_th = cat_th['pident']
                coverage_th = cat_th['coverage']
                evalue_min_th = cat_th['evalue_min']
                evalue_max_th = cat_th['evalue_max']
                print(f"     • Kategori özel: pident≥{pident_th:.1f}%, cov≥{coverage_th:.1f}%, "
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
                print(f"     • {dbtype}: {hits_after_filter} hit "
                      f"({avg_hits_per_intron:.2f} hit/intron ort.)")
        
        print(f"\n   📊 GENEL İSTATİSTİK:")
        print(f"     • Toplam intron: {total_introns_all}")
        print(f"     • Toplam hit: {total_hits_all}")
        print(f"     • Ortalama hit/intron: {total_hits_all/total_introns_all:.2f}")
        print(f"     • Maksimum hit/intron: 5 (kural uygulandı)")
        
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
            
            # Önce evalue range filtrele (MIN ve MAX arası)
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
            
            # Sonra coverage filtrele
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
            
            # TSV yaz
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
            
            # FASTA yaz
            out_fa = OUT_DIR / dbtype / "fasta" / f"{cat}.fa"
            out_fa.parent.mkdir(parents=True, exist_ok=True)
            with open(out_fa, "w") as f:
                for h in filtered_hits:
                    cid = h[0]
                    if cid in seqs:
                        for seq in seqs[cid][:1]:
                            f.write(f">{cid}|pident={h[2]}|evalue={h[10]}|bitscore={h[11]}\n{seq}\n")
            
            return filtered_hits
        except Exception as e:
            print(f"     ⚠️ BLAST hatası: {e}")
            return []
        finally:
            if tempfa.exists():
                tempfa.unlink()

# ===========================================
# 7. REPORT WRITER - GÜNCELLENDİ
# ===========================================
class ReportWriter:
    @staticmethod
    def write(thresholds, category_thresholds, cv_results, output_dir, summary_stats=None):
        print("\n[6/8] 📝 RAPOR HAZIRLANIYOR...")
        
        report_file = output_dir / "validation" / "optimization_report.txt"
        report_file.parent.mkdir(exist_ok=True)
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("BLAST OPTIMIZATION REPORT - ENHANCED v5.4\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("📁 INPUT FILES:\n")
            f.write(f"  Validator dir: {VALIDATOR_DIR}\n")
            f.write(f"  Database: {DB_PREMATURE.name}, {DB_MATURE.name}\n\n")
            
            f.write("🎯 OPTIMAL THRESHOLDS (VERİYE DAYALI):\n")
            f.write(f"  Percent Identity: ≥ {thresholds['pident']:.1f}%\n")
            f.write(f"  Coverage: ≥ {thresholds['coverage']:.1f}%\n")
            f.write(f"  Combined Score: ≥ {thresholds.get('combined_score', 0.7):.3f} (0.7*Identity + 0.3*Coverage)\n")
            f.write(f"  E-value MIN: ≥ {thresholds['evalue_min']:.1e} (10. percentile)\n")
            f.write(f"  E-value MAX: ≤ {thresholds['evalue_max']:.1e} (90. percentile)\n\n")
            
            f.write(f"📊 KORELASYON ANALİZİ:\n")
            f.write(f"  Identity-Coverage korelasyonu: {thresholds.get('correlation', 0):.3f}\n")
            if abs(thresholds.get('correlation', 0)) > 0.7:
                f.write(f"  ⚠️ YÜKSEK KORELASYON! Combined score kullanımı önerilir.\n")
            f.write("\n")
            
            if cv_results:
                f.write("🔬 CROSS-VALIDATION:\n")
                f.write(f"  pident mean: {cv_results['pident_mean']:.1f}%\n")
                f.write(f"  pident std: ±{cv_results['pident_std']:.2f}%\n")
                f.write(f"  stability score: {cv_results['stability_score']:.3f}\n\n")
            
            if category_thresholds:
                f.write("📊 CATEGORY-SPECIFIC THRESHOLDS:\n")
                for cat, th in category_thresholds.items():
                    f.write(f"  {cat}:\n")
                    f.write(f"    • pident ≥ {th['pident']:.1f}%\n")
                    f.write(f"    • coverage ≥ {th['coverage']:.1f}%\n")
                    f.write(f"    • combined score ≥ {th.get('combined_score', 0.7):.3f}\n")
                    f.write(f"    • evalue: {th['evalue_min']:.1e} - {th['evalue_max']:.1e} (5.-85. percentile)\n")
                    f.write(f"    • korelasyon: {th.get('correlation', 0):.3f}\n")
                f.write("\n")
            
            f.write("🎯 MULTIPLE MATCH FILTER:\n")
            f.write("  • Maximum hits per intron: 5\n")
            f.write("  • Selection criteria: Highest bitscore\n\n")
            
            if summary_stats:
                f.write("📊 FINAL BLAST STATISTICS:\n")
                f.write(f"  • Total introns: {summary_stats['total_introns']}\n")
                f.write(f"  • Total hits: {summary_stats['total_hits']}\n")
                f.write(f"  • Average hits per intron: {summary_stats['avg_hits_per_intron']:.2f}\n\n")
            
            f.write("⚙️  OPTIMIZATION PARAMETERS:\n")
            f.write("  • Exploratory BLAST: pident≥70, evalue≤0.01\n")
            f.write("  • Percentile threshold: 15th percentile (global)\n")
            f.write("  • Category threshold: 10th percentile\n")
            f.write("  • E-value MIN: 10th percentile (global), 5th percentile (category)\n")
            f.write("  • E-value MAX: 90th percentile (global), 85th percentile (category)\n")
            f.write("  • Combined score weights: 0.7 Identity + 0.3 Coverage\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("Report generated automatically by BLAST_DANEEL_ULTRA v5.4\n")
            f.write("=" * 80 + "\n")
        
        print(f"   ✅ Rapor kaydedildi: {report_file}")
        return report_file

# ===========================================
# 8. MANUAL BLAST
# ===========================================
class ManualBLAST:
    @staticmethod
    def run(pident, coverage, evalue_min, evalue_max):
        print("\n[1/2] 🎛️  MANUEL BLAST MODU")
        print("-" * 80)
        print(f"   • pident ≥ {pident}%")
        print(f"   • coverage ≥ {coverage}%")
        print(f"   • evalue: {evalue_min:.1e} - {evalue_max:.1e}")
        print(f"   • Maksimum hit/intron: 5")
        print("-" * 80)
        
        # Databaseleri kur
        ExploratoryBLAST._build_db(str(DB_PREMATURE), DB_DIR / "db_premature")
        ExploratoryBLAST._build_db(str(DB_MATURE), DB_DIR / "db_mature")
        
        thresholds = {
            'pident': pident,
            'coverage': coverage,
            'evalue_min': evalue_min,
            'evalue_max': evalue_max,
            'combined_score': 0.7  # Varsayılan
        }
        
        summary = FinalBLAST.run(thresholds, {})
        return summary

# ===========================================
# MAIN - GÜNCELLENDİ
# ===========================================
def main():
    print("\n" + "=" * 80)
    print("🚀 BLAST_DANEEL_ULTRA v5.4 - ENHANCED OPTIMIZATION!")
    print("=" * 80)
    
    # MODE 3: Manuel BLAST
    if mode == "3":
        summary = ManualBLAST.run(manual_pident, manual_coverage, manual_evalue_min, manual_evalue_max)
        
        print("\n📋 SUMMARY")
        print("-" * 80)
        with open(OUT_DIR / "summary.txt", "w") as f:
            for cat, db, hits, total_introns in summary:
                percentage = (hits / total_introns * 100) if total_introns > 0 else 0
                line = f"{cat}\t{db}\t{hits}\t{total_introns}\t{percentage:.1f}%"
                print(f"   {cat}\t{db}\t{hits} hits ({percentage:.1f}%)")
                f.write(line + "\n")
        
        print(f"\n✅ Sonuçlar: {OUT_DIR}/")
        return
    
    # MODE 1-2: Optimizasyon
    all_hits = ExploratoryBLAST.run()
    
    if len(all_hits) < 50:
        print("\n⚠️ Yeterli veri yok! Varsayılan eşikler kullanılacak.")
        thresholds = {
            'pident': 85.0,
            'coverage': 80.0,
            'combined_score': 0.7,
            'evalue_min': 1e-10,
            'evalue_max': 1e-2,
            'correlation': 0,
            'method': 'DEFAULT'
        }
        category_thresholds = {}
        cv_results = None
    else:
        # Feature extraction
        all_features, category_features = FeatureExtractor.extract(all_hits)
        
        # Global threshold optimization
        print("\n[3/8] 🧠 THRESHOLD OPTIMIZATION (VERİYE DAYALI + KOMBİNE SKOR)...")
        print("-" * 80)
        thresholds = ThresholdOptimizer.optimize(all_features)
        
        if not thresholds:
            thresholds = {
                'pident': 85.0,
                'coverage': 80.0,
                'combined_score': 0.7,
                'evalue_min': 1e-10,
                'evalue_max': 1e-2,
                'correlation': 0,
                'method': 'FALLBACK'
            }
        else:
            print(f"\n   ✅ OPTIMAL THRESHOLDS (VERİYE DAYALI):")
            print(f"      • pident ≥ {thresholds['pident']:.1f}%")
            print(f"      • coverage ≥ {thresholds['coverage']:.1f}%")
            print(f"      • combined score ≥ {thresholds['combined_score']:.3f}")
            print(f"      • evalue: {thresholds['evalue_min']:.1e} - {thresholds['evalue_max']:.1e}")
            print(f"      • {thresholds['method']}")
        
        # Category thresholds
        category_thresholds = {}
        for cat in CATS:
            if cat in category_features and len(category_features[cat]) >= 10:
                cat_th = ThresholdOptimizer.optimize_category(category_features[cat])
                if cat_th:
                    category_thresholds[cat] = cat_th
        
        # Cross-validation
        cv_results = CrossValidator.validate(all_features)
        
        # Plot
        PlotGenerator.generate(all_features, thresholds, OPT_DIR)
    
    # MODE 1: Final BLAST
    if mode == "1":
        summary = FinalBLAST.run(thresholds, category_thresholds)
        
        # Summary statistics
        total_introns = sum(item[3] for item in summary)
        total_hits = sum(item[2] for item in summary)
        
        summary_stats = {
            'total_introns': total_introns,
            'total_hits': total_hits,
            'avg_hits_per_intron': total_hits / total_introns if total_introns > 0 else 0
        }
        
        # Report
        ReportWriter.write(thresholds, category_thresholds, cv_results, OUT_DIR, summary_stats)
        
        print("\n[7/8] 📋 SUMMARY")
        print("-" * 80)
        
        with open(OUT_DIR / "summary.txt", "w") as f:
            for cat, db, hits, total_introns in summary:
                percentage = (hits / total_introns * 100) if total_introns > 0 else 0
                line = f"{cat}\t{db}\t{hits}\t{total_introns}\t{percentage:.1f}%"
                print(f"   {cat}\t{db}\t{hits} hits ({percentage:.1f}%)")
                f.write(line + "\n")
        
        print("\n" + "=" * 80)
        print("🎉 BLAST OPTIMIZATION COMPLETE! (ENHANCED v5.4)")
        print("=" * 80)
        print(f"\n📁 Output: {OUT_DIR}/")
        print(f"📊 Threshold plot: {OUT_DIR}/optimization/threshold_optimization.png")
        print(f"📑 Report: {OUT_DIR}/validation/optimization_report.txt")
        print(f"\n📌 E-value RANGE (VERİYE DAYALI): {thresholds['evalue_min']:.1e} - {thresholds['evalue_max']:.1e}")
        print(f"📌 COMBINED SCORE: {thresholds['combined_score']:.3f}")
        print("=" * 80)
    else:
        # MODE 2: Sadece optimizasyon
        ReportWriter.write(thresholds, category_thresholds, cv_results, OUT_DIR, None)
        
        print("\n" + "=" * 80)
        print("✅ OPTİMİZASYON TAMAMLANDI! (BLAST yapılmadı)")
        print(f"📁 Sonuçlar: {OUT_DIR}/")
        print("=" * 80)

if __name__ == "__main__":
    main()