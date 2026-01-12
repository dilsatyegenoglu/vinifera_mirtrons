#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
validator3_daneel_ultimate.py (REALTIME + FULL CLASS TSV + DEFAULT INPUT + ORIGINAL COLUMN ORDER FIXED + STATISTICS)
==================================================================================================================

Default input:
    rnafold_results/rna_inrange.tsv

Outputs:
    validator_results/
        validated_all.tsv
        validated_inrange.tsv
        invalid_entries.tsv
        validated_all.fa
        validated_inrange.fa
        classes/very_high.fa
        classes/high.fa
        classes/moderate.fa
        classes/low.fa
        classes/very_high.tsv
        classes/high.tsv
        classes/moderate.tsv
        classes/low.tsv
        summary.tsv
        statistics/ (new)
            gc_distribution.tsv
            mfe_distribution.tsv
            structural_features.tsv
            bps_analysis.tsv
            class_statistics.tsv

TSVs keep ORIGINAL column order + all new features appended at the END.
"""

import os
import csv
import time
import statistics
from collections import Counter, defaultdict

print("=== VALIDATOR3 DANEEL ULTIMATE (REALTIME FULL EDITION — ORIGINAL ORDER FIXED + STATISTICS) ===")

# -------------------------------------------------------------------
# DEFAULT INPUT
# -------------------------------------------------------------------
input_path = "rnafold_results/rna_inrange.tsv"
if not os.path.exists(input_path):
    raise FileNotFoundError(
        f"[ERROR] Cannot find: {input_path}\n"
        "Run RNAfold_DANEEL_PRO first!"
    )

print(f"[INFO] Using default input: {input_path}")

# -------------------------------------------------------------------
# LOAD TSV
# -------------------------------------------------------------------
def load_tsv(path):
    rows=[]
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

# -------------------------------------------------------------------
# BASIC CHECKS
# -------------------------------------------------------------------
def is_valid_seq(seq):
    return all(b in "ACGT" for b in seq)

def is_valid_dot(db):
    return all(c in "()." for c in db)

def gc_content(seq):
    return ((seq.count("G")+seq.count("C"))/len(seq))*100 if seq else 0

def basepair_count(db):
    return db.count("(")

def detect_bulges(db):
    bul=[]
    cur=0
    for c in db:
        if c==".":
            cur+=1
        else:
            if cur>0: bul.append(cur)
            cur=0
    if cur>0: bul.append(cur)
    return bul

def detect_stems(db):
    st=[]
    cur=0
    for c in db:
        if c=="(":
            cur+=1
        else:
            if cur>0: st.append(cur)
            cur=0
    if cur>0: st.append(cur)
    return st

def overhang_5(db):
    return len(db)-len(db.lstrip("."))

def overhang_3(db):
    return len(db)-len(db.rstrip("."))

# Plant BPS motifs
BPS_MOTIFS=["TACTAAC","CTAAC","CTRAC","CURA","TAC","AAC"]

def detect_branch(seq):
    seq=seq.upper()
    for m in BPS_MOTIFS:
        if m in seq:
            return True, m, seq.index(m)
    return False, "", -1

def amfe(mfe,length):
    return mfe/length if length else 0

def mfei(mfe,length,gc):
    if gc==0 or length==0:
        return 0
    return (mfe/length)/gc

# -------------------------------------------------------------------
# SCORE MODEL
# -------------------------------------------------------------------
def mirtron_score(r):
    score=0
    if 30<=r["GC"]<=65: score+=2
    if r["MFE"]<-25: score+=2
    if r["MFE"]<-35: score+=1
    if r["O5"]<=2: score+=2
    if r["O3"]<=2: score+=2
    if r["BPS"]: score+=2
    if r["MaxStem"]>=10: score+=2
    if r["BulgeCount"]<=3: score+=1
    return score

def score_class(s):
    if s>=10: return "very_high"
    if s>=7: return "high"
    if s>=5: return "moderate"
    if s>=3: return "low"
    return "reject"

# -------------------------------------------------------------------
# FASTA Writer
# -------------------------------------------------------------------
def write_fasta(path, rows):
    with open(path,"w") as f:
        for r in rows:
            hdr = f">{r['CoordID']}|len={r['Length']}|mfe={r['MFE']}|gc={r['GC']}|class={r['Class']}"
            f.write(hdr+"\n")
            seq=r["Sequence"]
            for i in range(0,len(seq),80):
                f.write(seq[i:i+80]+"\n")

# -------------------------------------------------------------------
# STATISTICS FUNCTIONS
# -------------------------------------------------------------------
def safe_float_convert(value, default=0.0):
    """Safely convert value to float, return default if conversion fails"""
    try:
        return float(value)
    except (ValueError, TypeError):
        return default

def safe_int_convert(value, default=0):
    """Safely convert value to int, return default if conversion fails"""
    try:
        return int(float(value))  # First convert to float to handle strings like "123.0"
    except (ValueError, TypeError):
        return default

def calculate_statistics(valid_rows):
    """Calculate comprehensive statistics for valid entries"""
    
    stats = {}
    
    # Convert string values to numeric for statistics
    gc_values = [safe_float_convert(r.get("GC", 0)) for r in valid_rows]
    mfe_values = [safe_float_convert(r.get("MFE", 0)) for r in valid_rows]
    length_values = [safe_int_convert(r.get("Length", 0)) for r in valid_rows]
    
    # GC Content Statistics
    stats["gc"] = {
        "mean": statistics.mean(gc_values) if gc_values else 0,
        "median": statistics.median(gc_values) if gc_values else 0,
        "stdev": statistics.stdev(gc_values) if len(gc_values) > 1 else 0,
        "min": min(gc_values) if gc_values else 0,
        "max": max(gc_values) if gc_values else 0,
        "distribution": Counter(round(gc) for gc in gc_values)
    }
    
    # MFE Statistics
    stats["mfe"] = {
        "mean": statistics.mean(mfe_values) if mfe_values else 0,
        "median": statistics.median(mfe_values) if mfe_values else 0,
        "stdev": statistics.stdev(mfe_values) if len(mfe_values) > 1 else 0,
        "min": min(mfe_values) if mfe_values else 0,
        "max": max(mfe_values) if mfe_values else 0,
        "distribution": Counter(round(mfe) for mfe in mfe_values)
    }
    
    # Length Statistics
    stats["length"] = {
        "mean": statistics.mean(length_values) if length_values else 0,
        "median": statistics.median(length_values) if length_values else 0,
        "min": min(length_values) if length_values else 0,
        "max": max(length_values) if length_values else 0
    }
    
    # Structural Features Statistics
    basepairs_values = [safe_int_convert(r.get("BasePairs", 0)) for r in valid_rows]
    bulge_count_values = [safe_int_convert(r.get("BulgeCount", 0)) for r in valid_rows]
    max_stem_values = [safe_int_convert(r.get("MaxStem", 0)) for r in valid_rows]
    overhang_5_values = [safe_int_convert(r.get("O5", 0)) for r in valid_rows]
    overhang_3_values = [safe_int_convert(r.get("O3", 0)) for r in valid_rows]
    
    stats["structural"] = {
        "basepairs_mean": statistics.mean(basepairs_values) if basepairs_values else 0,
        "bulge_count_mean": statistics.mean(bulge_count_values) if bulge_count_values else 0,
        "max_stem_mean": statistics.mean(max_stem_values) if max_stem_values else 0,
        "overhang_5_mean": statistics.mean(overhang_5_values) if overhang_5_values else 0,
        "overhang_3_mean": statistics.mean(overhang_3_values) if overhang_3_values else 0
    }
    
    # BPS Analysis
    bps_rows = [r for r in valid_rows if r.get("BPS") in [True, "True", "TRUE", "yes", "YES", "1", 1]]
    stats["bps"] = {
        "total_with_bps": len(bps_rows),
        "percentage": (len(bps_rows) / len(valid_rows)) * 100 if valid_rows else 0,
        "motif_frequency": Counter([r.get("BPS_motif", "unknown") for r in bps_rows]),
        "position_distribution": Counter([safe_int_convert(r.get("BPS_pos", -1)) for r in bps_rows if safe_int_convert(r.get("BPS_pos", -1)) != -1])
    }
    
    # Class-based Statistics
    class_stats = {}
    for class_name in ["very_high", "high", "moderate", "low", "reject"]:
        class_rows = [r for r in valid_rows if r.get("Class") == class_name]
        if class_rows:
            class_gc_values = [safe_float_convert(r.get("GC", 0)) for r in class_rows]
            class_mfe_values = [safe_float_convert(r.get("MFE", 0)) for r in class_rows]
            class_length_values = [safe_int_convert(r.get("Length", 0)) for r in class_rows]
            class_bps_count = len([r for r in class_rows if r.get("BPS") in [True, "True", "TRUE", "yes", "YES", "1", 1]])
            
            class_stats[class_name] = {
                "count": len(class_rows),
                "gc_mean": statistics.mean(class_gc_values) if class_gc_values else 0,
                "mfe_mean": statistics.mean(class_mfe_values) if class_mfe_values else 0,
                "length_mean": statistics.mean(class_length_values) if class_length_values else 0,
                "bps_percentage": (class_bps_count / len(class_rows)) * 100 if class_rows else 0
            }
    stats["classes"] = class_stats
    
    # Score Distribution
    score_values = [safe_float_convert(r.get("Score", 0)) for r in valid_rows]
    stats["scores"] = {
        "distribution": Counter(round(score) for score in score_values),
        "mean": statistics.mean(score_values) if score_values else 0,
        "median": statistics.median(score_values) if score_values else 0
    }
    
    return stats

def write_statistics(stats, output_dir, total_valid_entries):
    """Write statistics to TSV files"""
    stats_dir = os.path.join(output_dir, "statistics")
    os.makedirs(stats_dir, exist_ok=True)
    
    # GC Distribution
    with open(os.path.join(stats_dir, "gc_distribution.tsv"), "w") as f:
        f.write("GC_Value\tCount\n")
        for gc_val, count in sorted(stats["gc"]["distribution"].items()):
            f.write(f"{gc_val}\t{count}\n")
    
    # MFE Distribution
    with open(os.path.join(stats_dir, "mfe_distribution.tsv"), "w") as f:
        f.write("MFE_Value\tCount\n")
        for mfe_val, count in sorted(stats["mfe"]["distribution"].items()):
            f.write(f"{mfe_val}\t{count}\n")
    
    # Structural Features Summary
    with open(os.path.join(stats_dir, "structural_features.tsv"), "w") as f:
        f.write("Feature\tValue\n")
        f.write(f"Mean_BasePairs\t{stats['structural']['basepairs_mean']:.2f}\n")
        f.write(f"Mean_BulgeCount\t{stats['structural']['bulge_count_mean']:.2f}\n")
        f.write(f"Mean_MaxStem\t{stats['structural']['max_stem_mean']:.2f}\n")
        f.write(f"Mean_Overhang5\t{stats['structural']['overhang_5_mean']:.2f}\n")
        f.write(f"Mean_Overhang3\t{stats['structural']['overhang_3_mean']:.2f}\n")
    
    # BPS Analysis
    with open(os.path.join(stats_dir, "bps_analysis.tsv"), "w") as f:
        f.write("Metric\tValue\n")
        f.write(f"Total_With_BPS\t{stats['bps']['total_with_bps']}\n")
        f.write(f"Percentage_With_BPS\t{stats['bps']['percentage']:.2f}\n")
        f.write("\nBPS_Motif\tCount\n")
        for motif, count in stats['bps']['motif_frequency'].items():
            f.write(f"{motif}\t{count}\n")
    
    # Class Statistics
    with open(os.path.join(stats_dir, "class_statistics.tsv"), "w") as f:
        f.write("Class\tCount\tGC_Mean\tMFE_Mean\tLength_Mean\tBPS_Percentage\n")
        for class_name, class_data in stats["classes"].items():
            f.write(f"{class_name}\t{class_data['count']}\t{class_data['gc_mean']:.2f}\t{class_data['mfe_mean']:.2f}\t{class_data['length_mean']:.2f}\t{class_data['bps_percentage']:.2f}\n")
    
    # Overall Statistics Summary
    with open(os.path.join(stats_dir, "overall_statistics.tsv"), "w") as f:
        f.write("Metric\tValue\n")
        f.write(f"Total_Valid_Entries\t{total_valid_entries}\n")
        f.write(f"GC_Mean\t{stats['gc']['mean']:.2f}\n")
        f.write(f"GC_Median\t{stats['gc']['median']:.2f}\n")
        f.write(f"GC_StDev\t{stats['gc']['stdev']:.2f}\n")
        f.write(f"MFE_Mean\t{stats['mfe']['mean']:.2f}\n")
        f.write(f"MFE_Median\t{stats['mfe']['median']:.2f}\n")
        f.write(f"Length_Mean\t{stats['length']['mean']:.2f}\n")
        f.write(f"Length_Median\t{stats['length']['median']:.2f}\n")
        f.write(f"Score_Mean\t{stats['scores']['mean']:.2f}\n")
        f.write(f"Score_Median\t{stats['scores']['median']:.2f}\n")

# -------------------------------------------------------------------
# MAIN PIPELINE
# -------------------------------------------------------------------
rows = load_tsv(input_path)
total = len(rows)
print(f"[INFO] Loaded {total} introns\n")

valid=[]
invalid=[]
start = time.time()

for idx, r in enumerate(rows, start=1):

    cid=r["CoordID"]
    seq=r["Sequence"].strip()
    db=r["DotBracket"]
    
    # Safe conversion of MFE and Length
    try:
        mfe = float(r["MFE"])
    except (ValueError, TypeError):
        mfe = 0.0
        
    try:
        length = int(r["Length"])
    except (ValueError, TypeError):
        length = 0
        
    iso=r["Isoforms"]

    if not is_valid_seq(seq):
        r["Reason"]="Invalid sequence"
        invalid.append(r)
        continue

    if not is_valid_dot(db):
        r["Reason"]="Invalid dot-bracket"
        invalid.append(r)
        continue

    GC = gc_content(seq)
    BP = basepair_count(db)
    bul=detect_bulges(db)
    stems=detect_stems(db)
    O5=overhang_5(db)
    O3=overhang_3(db)
    bfound,bmotif,bpos = detect_branch(seq)
    AM=amfe(mfe,length)
    MI=mfei(mfe,length,GC)

    score = mirtron_score({
        "GC":GC,"MFE":mfe,"O5":O5,"O3":O3,
        "BPS":bfound,"MaxStem":max(stems) if stems else 0,
        "BulgeCount":len(bul)
    })
    c = score_class(score)

    # -------------------------------------------------------------------
    # ✔✔✔ ORIGINAL COLUMN ORDER FIX (CRITICAL)
    # -------------------------------------------------------------------
    rfeat = r.copy()     # ORIGINAL COLUMNS PRESERVED

    # New features appended at END
    rfeat.update({
        "GC": round(GC,3),
        "BasePairs": BP,
        "BulgeCount": len(bul),
        "MaxBulge": max(bul) if bul else 0,
        "MaxStem": max(stems) if stems else 0,
        "O5": O5,
        "O3": O3,
        "BPS": bfound,
        "BPS_motif": bmotif,
        "BPS_pos": bpos,
        "AMFE": round(AM,4),
        "MFEI": round(MI,4),
        "Score": score,
        "Class": c
    })

    valid.append(rfeat)

    print(f"[{idx:05d}/{total}] {cid} | MFE={mfe} | GC={round(GC,2)} | class={c}")

# -------------------------------------------------------------------
# OUTPUT STRUCTURE
# -------------------------------------------------------------------
outdir = "validator_results"
os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.join(outdir,"classes"), exist_ok=True)

def write_tsv(path, rows):
    if not rows: return
    with open(path,"w",newline="") as f:
        fn = list(rows[0].keys())    # ORIGINAL ORDER + appended features
        w  = csv.DictWriter(f,fieldnames=fn,delimiter="\t")
        w.writeheader()
        w.writerows(rows)

write_tsv(os.path.join(outdir,"invalid_entries.tsv"), invalid)
write_tsv(os.path.join(outdir,"validated_all.tsv"), valid)

valid_inrange = [r for r in valid if r["Class"]!="reject"]
write_tsv(os.path.join(outdir,"validated_inrange.tsv"), valid_inrange)

write_fasta(os.path.join(outdir,"validated_all.fa"), valid)
write_fasta(os.path.join(outdir,"validated_inrange.fa"), valid_inrange)

# CLASS FASTA & TSV
for cl in ["very_high","high","moderate","low"]:
    rows_cl = [r for r in valid if r["Class"] == cl]
    write_fasta(os.path.join(outdir,"classes",f"{cl}.fa"), rows_cl)
    write_tsv(os.path.join(outdir,"classes",f"{cl}.tsv"), rows_cl)

# STATISTICS
print("[INFO] Calculating comprehensive statistics...")
stats = calculate_statistics(valid)
write_statistics(stats, outdir, len(valid))  # total_valid_entries parametresi eklendi

# SUMMARY
summary={
    "Total": total,
    "Valid": len(valid),
    "Invalid": len(invalid),
    "VeryHigh": sum(r["Class"]=="very_high" for r in valid),
    "High": sum(r["Class"]=="high" for r in valid),
    "Moderate": sum(r["Class"]=="moderate" for r in valid),
    "Low": sum(r["Class"]=="low" for r in valid),
    "Reject": sum(r["Class"]=="reject" for r in valid)
}

with open(os.path.join(outdir,"summary.tsv"),"w") as f:
    f.write("Metric\tValue\n")
    for k,v in summary.items():
        f.write(f"{k}\t{v}\n")

elapsed = time.time() - start
print(f"\n=== VALIDATION COMPLETE in {elapsed:.1f} sec ===")
for k,v in summary.items():
    print(f"{k:12s}: {v}")

print("\nvalidator_results/ ready ✓")
print("Original column order preserved ✓")
print("Class FASTA + TSV generated ✓")
print("Comprehensive statistics generated ✓")
print("Pipeline-safe ✓")
print("REALTIME validator finished ✓")