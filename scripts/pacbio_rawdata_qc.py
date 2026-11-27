#!/usr/bin/env python3
"""
PacBio Revio Raw Data QC Report Generator
SMRT Link Run Details report 기준으로 HiFi raw data를 평가하여 HTML 리포트 생성

근거 문서: docs/103-604-700 SMRT Link Run Details report overview, Revio v13.3 Rev 01.pdf
적용 스펙 행: **15–20 kb 라이브러리** (문서 표에서 24-hr 와 짝지어진 행,
HiFi yield 100–120 Gb, Q30+ 90%). 실제 movie time 은 30 hr 이라 이 행보다
유리한 조건이다.

PacBio Revio v13.3 QC thresholds (15-20 kb WGS library 행):
  - HiFi Yield          : PASS ≥100 Gb | WARNING 70-100 Gb | FAIL <70 Gb
  - HiFi Read Length    : PASS ≥15 kb  | WARNING 10-15 kb  | FAIL <10 kb
  - HiFi Read Quality   : PASS ≥Q30   | WARNING Q20-Q30   | FAIL <Q20
  - Q30+ Bases          : PASS ≥90%   | WARNING 85-90%    | FAIL <85%
  - P1 (ZMW productive) : PASS 50-70% | WARNING 40-50% or 70-80% | FAIL <40% or >80%
  - Missing Adapter     : PASS ≤10%   | WARNING >10%      (문서가 주는 선은 10% 하나뿐)
  - Control Read Count  : PASS ≥500   | WARNING 250-499   | FAIL <250
  - Control Read Length : PASS ≥50 kb | WARNING 20-50 kb  | FAIL <20 kb
  - Control Concordance : PASS ≥88%   | WARNING <88%
  - Est. Coverage       : PASS ≥35×   | WARNING 25-35×    | FAIL <25×  (내부 기준)

control 세 지표와 P1 은 **장비가 계산해 둔 값을 읽는다**. 직접 계산하지 않는다 —
이유는 read_instrument_reports() 주석에 적었다.

Usage (single run):
  python3 pacbio_rawdata_qc.py --run-dir /data_4tb/pacbio_rawdata/r84285_XXXXXXXX \
      --sample-map 1_A01=TG-XXXX 1_B01=TG-XXXX 1_C01=TG-XXXX 1_D01=TG-XXXX \
      --output /path/to/output.html \
      [--movie-time 24]

Usage (multi-run, single project):
  python3 pacbio_rawdata_qc.py \
      --run-dir /data_4tb/pacbio_rawdata/r84285_XXXXXXXX /data_4tb/pacbio_rawdata/r84285_YYYYYYYY \
      --sample-map 1_A01=Sample1 1_B01=Sample2 1_D01=Sample3 2_A01=Sample4 \
      --output /path/to/output.html
"""

import argparse
import json
import os
import re
import glob
import sys
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime
from pathlib import Path


# ─── QC Thresholds (15-20 kb WGS 행) ────────────────────────────────────────
#
# `source` 는 그 숫자가 어디서 왔는지다. 섞지 않는다 — 제조사 문서가 말하지
# 않는 칸을 문서 기준인 척 두면, 나중에 "기준이 현실과 안 맞는다"는 잘못된
# 결론으로 이어진다 (Control RL 이 그 경로로 한 번 오해를 만들었다).
#   PACBIO  : 문서 103-604-700 (Revio v13.3) 에 적힌 값
#   INTERNAL: 우리 분석 요구에서 나온 값. 근거를 바로 밑에 적는다.
PACBIO = "PacBio 103-604-700 (Revio v13.3)"
INTERNAL = "TG 접합 검출 요구 — THRESHOLDS 주석 참고"

THRESHOLDS = {
    # 문서 스펙 표: 15-20 kb 라이브러리 -> 100-120 Gb / SMRT Cell.
    # 100 Gb 는 그 기대 구간의 아래끝이다.
    "hifi_yield_gb": {
        "pass_min": 100,
        "warn_min": 70,
        "label": "HiFi Yield (Gb)",
        "direction": "higher_better",
        "source": PACBIO,
    },
    # "Read lengths correspond to expected library insert size" -> 15-20 kb 행의 아래끝
    "hifi_read_length_mean_kb": {
        "pass_min": 15,
        "warn_min": 10,
        "label": "HiFi Read Length mean (kb)",
        "direction": "higher_better",
        "source": PACBIO,
    },
    "hifi_read_quality_median": {
        "pass_min": 30,
        "warn_min": 20,
        "label": "HiFi Read Quality median (QV)",
        "direction": "higher_better",
        "source": PACBIO,
    },
    # 스펙 표 15-20 kb 행의 Q30+ 값. (30 hr 행인 20-25 kb 는 85% 다 —
    # 우리는 15-20 kb 행을 쓰기로 했으므로 90% 를 유지한다.)
    "q30_base_pct": {
        "pass_min": 90,
        "warn_min": 85,
        "label": "Q30+ Bases (%)",
        "direction": "higher_better",
        "source": PACBIO,
    },
    # 문서: "a P1 in the range of ~50-70% is optimal", "<50% under-loaded",
    # "once P1 increases to >80% ... most likely entered an overloaded state"
    "p1_pct": {
        "pass_min": 50,
        "pass_max": 70,
        "warn_min": 40,
        "warn_max": 80,
        "label": "P1 – ZMW Productivity (%)",
        "direction": "range",   # 50-70% optimal
        "source": PACBIO,
    },
    # 문서가 말하는 것은 "Missing adapter should ideally be <10%" 한 줄뿐이다.
    # 예전에는 PASS <5% / WARNING 5-10% 로 쪼개 두었는데 5% 는 근거가 없었다.
    # 문서대로 10% 한 선만 쓰고, 그 아래를 우리가 다시 쪼개지 않는다.
    "missing_adapter_pct": {
        "pass_max": 10,
        "label": "Missing Adapter (%)",
        "direction": "lower_better_soft",
        "source": PACBIO,
    },
    "control_read_count": {
        "pass_min": 500,
        "warn_min": 250,
        "label": "Control Read Count",
        "direction": "higher_better",
        "source": PACBIO,
    },
    # "control read length is typically >=50 kb for 24 hrs movie time".
    # 우리는 30 hr 이므로 이 선은 보수적이다. 반드시 polymerase read length 로
    # 재야 한다 — read_instrument_reports() 주석 참고.
    "control_read_length_mean_kb": {
        "pass_min": 50,
        "warn_min": 20,
        "label": "Control Read Length mean (kb)",
        "direction": "higher_better",
        "source": PACBIO,
    },
    # 문서는 "mean & mode are typically ≥88%" 한 줄만 말한다. 그 아래를 다시
    # WARNING / FAIL 로 쪼갤 근거가 문서에 없으므로 쪼개지 않는다 (2단 판정).
    "control_concordance_pct": {
        "pass_min": 88,
        "label": "Control Read Concordance mean (%)",
        "direction": "higher_better_soft",
        "source": PACBIO,
    },
    # **제조사 문서는 이 기준을 줄 수 없다.** 문서의 "expected coverage 20x" 는
    # 인간 유전체(3.1 Gb) 기준의 estimate 이고, 우리는 마우스(2.7 Gb)에서
    # TG 통합 접합을 찾는다. 그래서 "TG 접합 검출에 몇 x 가 필요한가" 라는
    # 우리 분석 요구에서 나와야 한다.
    #
    # 근거:
    #   - 지금까지 납품한 분석은 35.3-39.3x 구간에서 수행됐다. 그 아래에서
    #     접합을 찾아낸 실적이 없다.
    #   - TG-6283 의 오른쪽 접합은 37.2x 에서 지지 리드가 ~1개였다. 즉 35x
    #     부근도 접합 하나를 리드 한 개로 버티는 지점이고, 더 낮아지면
    #     접합을 놓칠 위험이 실재한다.
    # 따라서 35x 를 "검증된 깊이", 25x 를 "접합 검출을 기대하기 어려운 깊이"로
    # 둔다. 이 숫자는 실적이 늘면 다시 볼 값이다.
    "est_coverage_x": {
        "pass_min": 35,
        "warn_min": 25,
        "label": "Est. Coverage (×)",
        "direction": "higher_better",
        "source": INTERNAL,
    },
}


# ─── Status evaluation ────────────────────────────────────────────────────────

def evaluate(metric_key, value):
    """Return 'PASS', 'WARNING', or 'FAIL' for a given metric value."""
    if value is None:
        return "N/A"
    t = THRESHOLDS[metric_key]
    direction = t["direction"]

    if direction == "higher_better":
        if value >= t["pass_min"]:
            return "PASS"
        elif value >= t["warn_min"]:
            return "WARNING"
        else:
            return "FAIL"

    elif direction == "higher_better_soft":
        # 문서가 기대값 한 개만 주는 지표. FAIL 선을 우리가 만들지 않는다.
        return "PASS" if value >= t["pass_min"] else "WARNING"

    elif direction == "lower_better_soft":
        return "PASS" if value <= t["pass_max"] else "WARNING"

    elif direction == "lower_better":
        if value <= t["pass_max"]:
            return "PASS"
        elif value <= t["warn_max"]:
            return "WARNING"
        else:
            return "FAIL"

    elif direction == "range":
        if t["pass_min"] <= value <= t["pass_max"]:
            return "PASS"
        elif t["warn_min"] <= value <= t["warn_max"]:
            return "WARNING"
        else:
            return "FAIL"

    return "N/A"


def overall_status(statuses):
    """Return worst status across a list of status strings."""
    if "FAIL" in statuses:
        return "FAIL"
    elif "WARNING" in statuses:
        return "WARNING"
    elif all(s == "PASS" for s in statuses if s != "N/A"):
        return "PASS"
    return "N/A"


# ─── Data parsing ─────────────────────────────────────────────────────────────

def parse_ccs_report(json_path):
    """Parse ccs_report.json and return a dict of metrics."""
    with open(json_path) as f:
        data = json.load(f)
    attrs = {a["id"]: a["value"] for a in data.get("attributes", [])}

    zmw_input       = attrs.get("ccs_processing.zmw_input", 0)
    zmw_passed      = attrs.get("ccs_processing.zmw_passed_yield", 0)
    lacking_passes  = attrs.get("ccs_processing.zmw_filtered_too_few_passes", 0)
    missing_adp     = attrs.get("ccs_processing.zmw_missing_adapters", 0)

    hifi_reads   = attrs.get("ccs_processing.number_of_ccs_reads_q20", 0)
    hifi_bases   = attrs.get("ccs_processing.total_number_of_ccs_bases_q20", 0)
    hifi_len_mean= attrs.get("ccs_processing.mean_ccs_readlength_q20", 0)
    hifi_len_n50 = attrs.get("ccs_processing.ccs_readlength_n50_q20", 0)
    hifi_qv_med  = attrs.get("ccs_processing.median_qv_q20", 0)
    q30_pct      = attrs.get("ccs_processing.base_percentage_q30", 0)
    mean_passes  = attrs.get("ccs_processing.mean_npasses_q20", 0)
    lq_reads     = attrs.get("ccs_processing.number_of_ccs_reads_lq", 0)

    # P1 proxy: ZMWs that had productive signal (≈ zmw_input - P0 empty ZMWs).
    # "lacking full passes" ≈ P0 (empty/unloaded ZMWs with no usable signal).
    # This is an approximation; exact P1 requires SMRT Link instrument metrics.
    p1_zmws = zmw_input - lacking_passes
    p1_pct = (p1_zmws / zmw_input * 100) if zmw_input > 0 else None
    missing_pct = (missing_adp / zmw_input * 100) if zmw_input > 0 else None

    return {
        "zmw_input": zmw_input,
        "zmw_passed": zmw_passed,
        "hifi_reads": hifi_reads,
        "lq_reads": lq_reads,
        "hifi_yield_gb": hifi_bases / 1e9,
        "hifi_read_length_mean_kb": hifi_len_mean / 1e3,
        "hifi_read_length_n50_kb": hifi_len_n50 / 1e3,
        "hifi_read_quality_median": hifi_qv_med,
        "q30_base_pct": q30_pct,
        "mean_passes": mean_passes,
        "p1_pct": p1_pct,
        "missing_adapter_pct": missing_pct,
    }


def read_instrument_reports(stats_dir, movie_name):
    """장비가 이미 계산해 둔 control / loading 지표를 읽는다.

    **왜 직접 계산하지 않나.** 예전에는 control 지표를
    `metadata/*.sequencing_control.subreads.bam` 을 열어 SEQ 길이 평균으로
    계산했다. 그 값은 control 템플릿을 **한 번 통과한 subread 의 길이**
    (23~25 kb) 다. 문서 103-604-700 이 말하는 `control read length is
    typically ≥50 kb` 는 control 의 **polymerase read length** 다 — 한 ZMW 의
    모든 통과를 합친 길이. 단위가 달라서 28개 샘플 전부가 WARNING 으로
    찍혔고, 임계값이 틀린 것으로 오해했다. 실제로는 임계값이 맞고 측정이
    틀렸다 (장비값으로 재집계하면 과거 24웰 중 23웰 PASS).

    read count 도 같은 문제였다 — subread 개수(수천)를 세고 있었고 장비가
    말하는 control read 수는 ZMW 기준이다.

    P1 도 마찬가지다. ccs_report 기반 계산은 근사치이고, 원래 주석도
    "exact P1 requires SMRT Link instrument metrics" 라고 적어 두었다.
    그 instrument metrics 가 같은 statistics/ 폴더 안에 있었다.

    읽는 곳: statistics/<movie>.reports.zip 안의
        control.report.json  -> reads_n / readlength_mean / concordance_mean·mode
        loading.report.json  -> productivity_0_n / _1_n / _2_n

    reports.zip 이 없으면 전부 None 을 돌려준다. subreads BAM 으로 추정해서
    틀린 단위의 값을 임계값에 들이대지 않는다.
    """
    out = {
        "control_read_count": None,
        "control_read_length_mean_kb": None,
        "control_concordance_pct": None,
        "control_concordance_mode_pct": None,
        "p1_pct_instrument": None,
        "p2_pct_instrument": None,
    }
    zip_path = os.path.join(stats_dir, f"{movie_name}.reports.zip")
    if not os.path.exists(zip_path):
        return out

    try:
        with zipfile.ZipFile(zip_path) as z:
            names = set(z.namelist())

            if "control.report.json" in names:
                attrs = {
                    a["id"]: a["value"]
                    for a in json.loads(z.read("control.report.json")).get("attributes", [])
                }
                rl = attrs.get("control.readlength_mean")
                conc = attrs.get("control.concordance_mean")
                mode = attrs.get("control.concordance_mode")
                out["control_read_count"] = attrs.get("control.reads_n")
                out["control_read_length_mean_kb"] = rl / 1e3 if rl is not None else None
                out["control_concordance_pct"] = conc * 100 if conc is not None else None
                out["control_concordance_mode_pct"] = mode * 100 if mode is not None else None

            if "loading.report.json" in names:
                attrs = {
                    a["id"].rsplit(".", 1)[-1]: a["value"]
                    for a in json.loads(z.read("loading.report.json")).get("attributes", [])
                }
                p0 = attrs.get("productivity_0_n")
                p1 = attrs.get("productivity_1_n")
                p2 = attrs.get("productivity_2_n")
                if None not in (p0, p1, p2):
                    total = p0 + p1 + p2
                    if total > 0:
                        out["p1_pct_instrument"] = p1 / total * 100
                        out["p2_pct_instrument"] = p2 / total * 100
    except Exception as e:
        print(f"  [WARN] {zip_path} 를 읽을 수 없다: {e}", file=sys.stderr)

    return out


def get_barcode(hifi_reads_dir):
    """Return barcode string from hifi_reads BAM filename (optional dir)."""
    if not os.path.isdir(hifi_reads_dir):
        return "unknown"
    bams = glob.glob(os.path.join(hifi_reads_dir, "*.hifi_reads.bc*.bam"))
    bams = [b for b in bams if not b.endswith(".pbi")]
    if bams:
        m = re.search(r"(bc\d+)", os.path.basename(bams[0]))
        return m.group(1) if m else "unknown"
    return "unknown"


def read_declared_insert(well_dir, movie_name):
    """랩이 run 셋업 때 입력한 라이브러리 insert 크기를 읽는다.

    <well>/metadata/<movie>.metadata.xml 의 <InsertSize> — sizing 측정 후
    사람이 넣는 값이다 (반올림되지 않은 실측값이 들어온다. 예: 18519).

    **왜 리포트에 넣나.** 실측 HiFi read length 는 이 선언값보다 **항상**
    작다 — 28웰 전수 집계에서 예외가 0건이었고, 비율은 62.4-83.6%
    (중위수 74.0%, 표준편차 5.3%p) 였다. 평균뿐 아니라 N50 도 전부 선언값
    아래이므로 분포가 통째로 내려앉는 것이고 꼬리 효과가 아니다.

    그래서 Read Length 가 낮게 나왔을 때 원인이 두 갈래다 —
    목표 insert 를 낮게 잡았나, 아니면 라이브러리가 선언 대비 못 나왔나.
    선언값을 같이 적어 두지 않으면 나중에 이 둘을 구분할 수 없다.
    랩이 run 마다 의도적으로 바꾸는 값이라 더욱 기록이 필요하다
    (2026-11 run 은 이 차이를 감안해 18.1 -> 21.1 kb 로 올렸다).

    비율에는 **판정을 붙이지 않는다.** run 안에서는 2-4%p 로 모이고 run
    사이에서 72% vs 82% 로 갈리는, 샘플이 아니라 prep 배치의 성질이다.
    선을 그으면 배치 차이를 샘플 경고로 잘못 보고하게 된다.
    """
    out = {"declared_insert_bp": None, "library_type": None, "size_selection": None}
    xml_path = os.path.join(well_dir, "metadata", f"{movie_name}.metadata.xml")
    if not os.path.exists(xml_path):
        return out
    try:
        root = ET.parse(xml_path).getroot()
    except ET.ParseError as e:
        print(f"  [WARN] {xml_path} 를 읽을 수 없다: {e}", file=sys.stderr)
        return out

    def first_text(name):
        for el in root.iter():
            if el.tag.rsplit("}", 1)[-1] == name and (el.text or "").strip():
                return el.text.strip()
        return None

    raw = first_text("InsertSize")
    if raw is not None:
        try:
            out["declared_insert_bp"] = int(float(raw))
        except ValueError:
            print(f"  [WARN] InsertSize 를 숫자로 읽을 수 없다: {raw!r}", file=sys.stderr)
    out["library_type"] = first_text("LibraryType")
    out["size_selection"] = first_text("SizeSelectionEnabled")
    return out


def collect_well_metrics(run_dir, well, sample_name, genome_size_gb=0.0):
    """Collect all QC metrics for one well.

    genome_size_gb 를 주면 Est. Coverage 를 계산하고 판정까지 한다 (0 이면
    coverage 는 N/A). coverage 는 유전체 크기를 알아야 나오는 값이라서
    여기서만 계산할 수 있다.
    """
    well_dir = os.path.join(run_dir, well)
    stats_dir = os.path.join(well_dir, "statistics")

    ccs_jsons = glob.glob(os.path.join(stats_dir, "*.ccs_report.json"))
    if not ccs_jsons:
        print(f"  [WARN] No ccs_report.json found for {well}", file=sys.stderr)
        return None

    movie_name = os.path.basename(ccs_jsons[0]).replace(".ccs_report.json", "")
    metrics = parse_ccs_report(ccs_jsons[0])
    metrics["movie"] = movie_name
    metrics["well"] = well
    metrics["run"] = os.path.basename(os.path.normpath(run_dir))
    metrics["sample"] = sample_name
    metrics["barcode"] = get_barcode(os.path.join(well_dir, "hifi_reads"))

    inst = read_instrument_reports(stats_dir, movie_name)
    metrics.update({
        "control_read_count": inst["control_read_count"],
        "control_read_length_mean_kb": inst["control_read_length_mean_kb"],
        "control_concordance_pct": inst["control_concordance_pct"],
        "control_concordance_mode_pct": inst["control_concordance_mode_pct"],
    })

    # P1: 장비 실측값이 있으면 그것을 쓰고, ccs_report 근사치는 비교용으로 남긴다.
    metrics["p1_pct_ccs_proxy"] = metrics.get("p1_pct")
    if inst["p1_pct_instrument"] is not None:
        metrics["p1_pct"] = inst["p1_pct_instrument"]
        metrics["p1_source"] = "instrument"
    else:
        metrics["p1_source"] = "ccs_proxy"
    metrics["p2_pct"] = inst["p2_pct_instrument"]

    if inst["control_read_length_mean_kb"] is None:
        print(f"    [INFO] {well}: {movie_name}.reports.zip 에서 control 지표를 "
              f"읽지 못했다 — N/A 로 둔다 (subreads BAM 으로 추정하지 않는다)",
              file=sys.stderr)

    # 선언 insert 크기와 실측 read length 의 비. 판정은 붙이지 않는다 (진단용).
    metrics.update(read_declared_insert(well_dir, movie_name))
    di = metrics.get("declared_insert_bp")
    rl = metrics.get("hifi_read_length_mean_kb")
    metrics["rl_insert_pct"] = (
        rl * 1e3 / di * 100 if (di and di > 0 and rl is not None) else None
    )

    # Est. Coverage. 유전체 크기를 모르면 판정하지 않는다.
    if genome_size_gb and genome_size_gb > 0 and metrics.get("hifi_yield_gb") is not None:
        metrics["est_coverage_x"] = metrics["hifi_yield_gb"] / genome_size_gb
    else:
        metrics["est_coverage_x"] = None
    metrics["genome_size_gb"] = genome_size_gb

    # evaluate statuses
    status_keys = [
        "hifi_yield_gb", "est_coverage_x",
        "hifi_read_length_mean_kb", "hifi_read_quality_median",
        "q30_base_pct", "p1_pct", "missing_adapter_pct",
        "control_read_count", "control_read_length_mean_kb",
        "control_concordance_pct",
    ]
    metrics["statuses"] = {k: evaluate(k, metrics.get(k)) for k in status_keys}
    metrics["overall"] = overall_status(list(metrics["statuses"].values()))

    return metrics


# ─── HTML generation ──────────────────────────────────────────────────────────

STATUS_STYLE = {
    "PASS":    ("✅", "#d4edda", "#155724"),
    "WARNING": ("⚠️",  "#fff3cd", "#856404"),
    "FAIL":    ("❌", "#f8d7da", "#721c24"),
    "N/A":     ("—",  "#f8f9fa", "#6c757d"),
}


def badge(status):
    icon, bg, fg = STATUS_STYLE.get(status, ("?", "#eee", "#333"))
    return (
        f'<span style="background:{bg};color:{fg};padding:2px 8px;'
        f'border-radius:4px;font-weight:bold;white-space:nowrap">'
        f'{icon} {status}</span>'
    )


def fmt(value, decimals=1, suffix=""):
    if value is None:
        return "N/A"
    return f"{value:.{decimals}f}{suffix}"


def pct(value, decimals=1):
    """None 이면 'N/A' — 'N/A%' 가 되지 않게 한다."""
    return "N/A" if value is None else f"{value:.{decimals}f}%"


def build_summary_table(all_metrics, genome_size_gb):
    def est_coverage(m):
        if m.get("est_coverage_x") is None:
            return "N/A"
        return f"{m['est_coverage_x']:.1f}×"

    cols = [
        ("Sample",                lambda m: f"<b>{m['sample']}</b>"),
        ("Run",                   lambda m: f'<small>{m["run"]}</small>'),
        ("Well",                  lambda m: m["well"]),
        ("Barcode",               lambda m: m["barcode"]),
        ("Overall",               lambda m: badge(m["overall"])),
        ("HiFi Reads (M)",        lambda m: fmt(m["hifi_reads"] / 1e6, 2)),
        ("HiFi Yield (Gb)",       lambda m: f"{fmt(m['hifi_yield_gb'])} {badge(m['statuses']['hifi_yield_gb'])}"),
        (f"Est. Coverage<br><small>({genome_size_gb} Gb genome · 내부 기준)</small>",
                                  lambda m: f"<b>{est_coverage(m)}</b> {badge(m['statuses']['est_coverage_x'])}"),
        ("Declared Insert (kb)<br><small>(랩 입력값 · 판정 없음)</small>",
                                  lambda m: ("N/A" if m.get("declared_insert_bp") is None
                                             else f"{m['declared_insert_bp'] / 1e3:.1f}")),
        ("Read Length mean (kb)", lambda m: f"{fmt(m['hifi_read_length_mean_kb'])} {badge(m['statuses']['hifi_read_length_mean_kb'])}"),
        ("Read Length N50 (kb)",  lambda m: fmt(m["hifi_read_length_n50_kb"])),
        ("RL / Insert<br><small>(진단용 · 판정 없음)</small>",
                                  lambda m: pct(m.get("rl_insert_pct"))),
        ("Read Quality (median)", lambda m: f"Q{m['hifi_read_quality_median']} {badge(m['statuses']['hifi_read_quality_median'])}"),
        ("Q30+ Bases (%)",        lambda m: f"{fmt(m['q30_base_pct'])}% {badge(m['statuses']['q30_base_pct'])}"),
        ("P1 (%)",                lambda m: f"{fmt(m['p1_pct'])}% {badge(m['statuses']['p1_pct'])}"),
        ("P2 (%)",                lambda m: pct(m.get("p2_pct"), 2)),
        ("Missing Adapter (%)",   lambda m: f"{fmt(m['missing_adapter_pct'])}% {badge(m['statuses']['missing_adapter_pct'])}"),
        ("Mean Passes",           lambda m: fmt(m["mean_passes"], 1)),
        ("Control Reads",         lambda m: f"{int(m['control_read_count']) if m['control_read_count'] is not None else 'N/A'} {badge(m['statuses']['control_read_count'])}"),
        ("Control RL mean (kb)",  lambda m: f"{fmt(m['control_read_length_mean_kb'])} {badge(m['statuses']['control_read_length_mean_kb'])}"),
        ("Control Concordance<br><small>(mean / mode)</small>",
                                  lambda m: f"{pct(m.get('control_concordance_pct'))} / "
                                            f"{pct(m.get('control_concordance_mode_pct'))} "
                                            f"{badge(m['statuses']['control_concordance_pct'])}"),
    ]

    header = "".join(f"<th>{c[0]}</th>" for c in cols)
    rows = []
    for m in all_metrics:
        cells = "".join(f"<td>{fn(m)}</td>" for _, fn in cols)
        bg = STATUS_STYLE.get(m["overall"], ("", "#fff", ""))[1]
        rows.append(f'<tr style="background:{bg}">{cells}</tr>')

    return f"""
<div class="table-scroll">
<table class="qc-table">
  <thead><tr>{header}</tr></thead>
  <tbody>{"".join(rows)}</tbody>
</table>
</div>"""


def build_threshold_table():
    rows = []
    for key, t in THRESHOLDS.items():
        d = t["direction"]
        if d == "higher_better":
            crit = f"PASS ≥{t['pass_min']} | WARNING ≥{t['warn_min']} | FAIL &lt;{t['warn_min']}"
        elif d == "lower_better":
            crit = f"PASS ≤{t['pass_max']}% | WARNING ≤{t['warn_max']}% | FAIL >{t['warn_max']}%"
        elif d == "higher_better_soft":
            crit = f"PASS ≥{t['pass_min']}% | WARNING &lt;{t['pass_min']}%"
        elif d == "lower_better_soft":
            crit = f"PASS ≤{t['pass_max']}% | WARNING >{t['pass_max']}%"
        elif d == "range":
            crit = (f"PASS {t['pass_min']}–{t['pass_max']}% | "
                    f"WARNING {t['warn_min']}–{t['pass_min']}% or {t['pass_max']}–{t['warn_max']}% | "
                    f"FAIL &lt;{t['warn_min']}% or >{t['warn_max']}%")
        src = t.get("source", "—")
        tag = ('<span style="background:#e7f1ff;color:#0b5394;padding:1px 6px;'
               'border-radius:3px;font-size:0.85em;white-space:nowrap">문서</span>'
               if src == PACBIO else
               '<span style="background:#fde9d9;color:#8a4b08;padding:1px 6px;'
               'border-radius:3px;font-size:0.85em;white-space:nowrap">내부</span>')
        rows.append(f"<tr><td><b>{t['label']}</b></td><td>{crit}</td>"
                    f"<td>{tag} <small>{src}</small></td></tr>")
    return f"""
<table class="qc-table">
  <thead><tr><th>Metric</th><th>QC Criteria (15–20 kb WGS 행, Revio v13.3)</th>
  <th>기준 출처</th></tr></thead>
  <tbody>{"".join(rows)}</tbody>
</table>
<p style="font-size:0.9em;color:#555;margin-top:6px">
  <b>문서</b> = PacBio 103-604-700 (SMRT Link Run Details report overview, Revio v13.3) 에
  적힌 값. <b>내부</b> = 제조사 문서가 줄 수 없어 우리 분석 요구에서 정한 값 —
  근거는 스크립트의 THRESHOLDS 주석에 적어 두었다.
  문서가 기대값 하나만 말하는 지표(Missing Adapter, Control Concordance)는
  그 아래를 임의로 쪼개지 않고 2단(PASS/WARNING)으로만 판정한다.
</p>"""


def build_bar_chart_js(all_metrics):
    """Build inline Chart.js code for key metrics bar charts."""
    samples = [m["sample"] for m in all_metrics]
    colors = {
        "PASS":    "rgba(40,167,69,0.7)",
        "WARNING": "rgba(255,193,7,0.7)",
        "FAIL":    "rgba(220,53,69,0.7)",
        "N/A":     "rgba(108,117,125,0.5)",
    }

    def bar_dataset(label, values, statuses):
        bg = [colors.get(s, colors["N/A"]) for s in statuses]
        return json.dumps({"label": label, "data": values, "backgroundColor": bg})

    charts = []
    chart_configs = [
        ("hifi_yield_gb",             "HiFi Yield (Gb)",          "hifi_yield_gb"),
        ("hifi_read_length_mean_kb",  "HiFi Read Length mean (kb)", "hifi_read_length_mean_kb"),
        ("hifi_read_quality_median",  "HiFi Read Quality (QV)",   "hifi_read_quality_median"),
        ("q30_base_pct",              "Q30+ Bases (%)",            "q30_base_pct"),
        ("p1_pct",                    "P1 – ZMW Productivity (%)", "p1_pct"),
        ("control_read_count",        "Control Read Count",        "control_read_count"),
    ]

    for i, (key, label, status_key) in enumerate(chart_configs):
        vals = [round(m.get(key) or 0, 2) for m in all_metrics]
        sts  = [m["statuses"].get(status_key, "N/A") for m in all_metrics]
        ds   = bar_dataset(label, vals, sts)
        charts.append(f"""
  new Chart(document.getElementById('chart_{i}'), {{
    type: 'bar',
    data: {{
      labels: {json.dumps(samples)},
      datasets: [{ds}]
    }},
    options: {{
      responsive: true,
      plugins: {{
        legend: {{ display: false }},
        title: {{ display: true, text: '{label}', font: {{ size: 14 }} }}
      }},
      scales: {{ y: {{ beginAtZero: true }} }}
    }}
  }});""")

    canvas_html = "".join(
        f'<div class="chart-box"><canvas id="chart_{i}"></canvas></div>'
        for i in range(len(chart_configs))
    )
    js = "<script>\n" + "\n".join(charts) + "\n</script>"
    return canvas_html, js


def build_html(all_metrics, run_dirs, movie_time, genome_size_gb):
    if isinstance(run_dirs, list):
        run_names = [os.path.basename(os.path.normpath(d)) for d in run_dirs]
    else:
        run_names = [os.path.basename(os.path.normpath(run_dirs))]
    run_label = " + ".join(run_names)
    title_run = run_names[0] if len(run_names) == 1 else f"{len(run_names)} runs"
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    summary_table = build_summary_table(all_metrics, genome_size_gb)
    threshold_table = build_threshold_table()
    canvas_html, charts_js = build_bar_chart_js(all_metrics)

    pass_ct    = sum(1 for m in all_metrics if m["overall"] == "PASS")
    warn_ct    = sum(1 for m in all_metrics if m["overall"] == "WARNING")
    fail_ct    = sum(1 for m in all_metrics if m["overall"] == "FAIL")

    return f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<title>PacBio HiFi Raw QC – {title_run}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4/dist/chart.umd.min.js"></script>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; background: #f5f7fa; color: #212529; }}
  .header {{ background: linear-gradient(135deg,#1a1a2e,#16213e); color: white; padding: 24px 32px; }}
  .header h1 {{ margin: 0 0 4px; font-size: 1.6em; }}
  .header p  {{ margin: 0; opacity: 0.75; font-size: 0.9em; }}
  .content {{ padding: 24px 32px; max-width: 1400px; margin: auto; }}
  h2 {{ color: #1a1a2e; border-bottom: 2px solid #dee2e6; padding-bottom: 6px; margin-top: 32px; }}
  .kpi-row {{ display: flex; gap: 16px; margin: 20px 0; flex-wrap: wrap; }}
  .kpi {{ background: white; border-radius: 8px; padding: 16px 24px; flex: 1; min-width: 120px;
           box-shadow: 0 1px 4px rgba(0,0,0,.08); text-align: center; }}
  .kpi .num {{ font-size: 2em; font-weight: bold; }}
  .kpi .lbl {{ font-size: 0.85em; color: #6c757d; margin-top: 4px; }}
  .kpi.pass  .num {{ color: #155724; }}
  .kpi.warn  .num {{ color: #856404; }}
  .kpi.fail  .num {{ color: #721c24; }}
  .qc-table {{ border-collapse: collapse; width: 100%; background: white;
               border-radius: 8px; overflow: hidden;
               box-shadow: 0 1px 4px rgba(0,0,0,.08); font-size: 0.88em; }}
  .qc-table th {{ background: #1a1a2e; color: white; padding: 10px 12px; text-align: left; white-space: nowrap; }}
  .qc-table td {{ padding: 9px 12px; border-bottom: 1px solid #dee2e6; vertical-align: middle; }}
  .qc-table tr:last-child td {{ border-bottom: none; }}
  .charts-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 20px; margin-top: 16px; }}
  .chart-box {{ background: white; border-radius: 8px; padding: 16px;
                box-shadow: 0 1px 4px rgba(0,0,0,.08); }}
  .table-scroll {{ overflow-x: auto; -webkit-overflow-scrolling: touch;
                   border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,.08); }}
  .table-scroll .qc-table {{ box-shadow: none; border-radius: 0; }}
  .note {{ background: #e8f4f8; border-left: 4px solid #0066cc; padding: 12px 16px;
           border-radius: 4px; font-size: 0.88em; margin: 16px 0; }}
  .ts-section {{ margin-bottom: 28px; }}
  .ts-section h3 {{ color: #1a1a2e; font-size: 1.05em; margin: 0 0 10px;
                    display: flex; align-items: center; gap: 8px; }}
  .ts-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(420px, 1fr)); gap: 14px; }}
  .ts-card {{ background: white; border-radius: 8px; padding: 16px 18px;
              box-shadow: 0 1px 4px rgba(0,0,0,.08); border-left: 4px solid #dee2e6; }}
  .ts-card.warn {{ border-left-color: #ffc107; }}
  .ts-card.fail {{ border-left-color: #dc3545; }}
  .ts-card.info {{ border-left-color: #0066cc; }}
  .ts-card h4 {{ margin: 0 0 8px; font-size: 0.95em; color: #333; }}
  .ts-card ul {{ margin: 0; padding-left: 18px; font-size: 0.88em; color: #444; line-height: 1.7; }}
  .ts-card .cause {{ font-size: 0.82em; color: #888; margin-top: 6px; font-style: italic; }}
  .flowchart {{ background: white; border-radius: 8px; padding: 20px 24px;
                box-shadow: 0 1px 4px rgba(0,0,0,.08); font-size: 0.88em; }}
  .flow-step {{ display: flex; align-items: flex-start; gap: 14px; margin-bottom: 14px; }}
  .flow-num {{ background: #1a1a2e; color: white; border-radius: 50%; width: 26px; height: 26px;
               display: flex; align-items: center; justify-content: center;
               font-weight: bold; flex-shrink: 0; font-size: 0.85em; }}
  .flow-body {{ flex: 1; }}
  .flow-body b {{ color: #1a1a2e; }}
  footer {{ text-align: center; padding: 24px; color: #aaa; font-size: 0.8em; }}
</style>
</head>
<body>
<div class="header">
  <h1>PacBio HiFi Raw Data QC Report</h1>
  <p>Run: {run_label} &nbsp;|&nbsp; Movie time: <b>{movie_time} hr</b> (실제)
     &nbsp;|&nbsp; Generated: {now}</p>
</div>
<div class="content">

<h2>Overall Summary</h2>
<div class="kpi-row">
  <div class="kpi"><div class="num">{len(all_metrics)}</div><div class="lbl">Total Samples</div></div>
  <div class="kpi pass"><div class="num">{pass_ct}</div><div class="lbl">PASS</div></div>
  <div class="kpi warn"><div class="num">{warn_ct}</div><div class="lbl">WARNING</div></div>
  <div class="kpi fail"><div class="num">{fail_ct}</div><div class="lbl">FAIL</div></div>
</div>

<div class="note">
  <b>기준:</b> PacBio Revio v13.3 스펙 표의 <b>15–20 kb 라이브러리 행</b>
  (HiFi yield 100–120 Gb, Q30+ 90%) — SMRT Link Run Details report overview,
  PN 103-604-700 Rev 01.<br>
  <b>movie time:</b> 이 run 의 실제 collection time 은 <b>{movie_time} hr</b> 입니다.
  문서 표는 15–20 kb 행을 24 hr 과 짝지었고 30 hr 은 20–25 kb 행에 둡니다 —
  즉 <b>우리 조건과 정확히 일치하는 행이 문서에 없습니다.</b> 더 긴 movie 는
  yield 와 read length 에 유리하므로 15–20 kb 행을 쓰는 것은 보수적인 선택입니다
  (Control Read Length ≥50 kb 도 24 hr 기준이라 같은 의미로 보수적입니다).<br>
  <b>Est. Coverage</b> = HiFi Yield / {genome_size_gb} Gb. 이 지표의 임계값만은
  제조사 문서가 아니라 <b>내부 기준</b>입니다 — 아래 기준표의 출처 열을 보십시오.
</div>

<div class="note">
  <b>Declared Insert 과 RL / Insert 읽는 법 (판정 없는 진단 지표):</b><br>
  Declared Insert 는 랩이 sizing 후 run 셋업에 입력한 값이고, 실측 HiFi read
  length 는 <b>이 값보다 항상 작습니다</b> — 자체 6 run 28웰 전수 집계에서
  예외가 0건이었습니다. 평균뿐 아니라 N50 도 28/28 선언값 아래여서, 짧은 조각이
  평균을 끌어내린 게 아니라 분포가 통째로 내려앉습니다.
  <ul style="margin:6px 0">
    <li>관측 비율 <b>62.4 – 83.6%</b> (중위수 74.0%, 표준편차 5.3%p)</li>
    <li>run 안에서는 2–4%p 로 모이고 <b>run 사이에서 72% 대 82% 로 갈립니다</b> —
        샘플 DNA 품질보다 prep 배치/시약의 성질입니다. 그래서 이 비율에는
        합불 판정을 붙이지 않았습니다.</li>
    <li>선언 크기 자체는 비율과 거의 무관합니다 (r = +0.17) — 크게 sizing 해도
        깎이는 비율이 더 커지지는 않습니다.</li>
  </ul>
  쓰는 법: Read Length 가 ⚠️ 인데 이 비율이 정상(≈75% 이상)이면 라이브러리가
  아니라 <b>목표 insert 가 낮았던</b> 것입니다. 비율까지 낮으면 prep 을 봐야
  합니다. 평균 15 kb 를 받으려면 선언 insert 는
  <b>최악 run 비율(71.4%)에서도 통과하려면 22 kb</b> 가 필요합니다
  (74.8% 평균 기준으로는 20 kb, 이번 run 은 21.1 kb 를 선언해 4/4 통과).
</div>

<h2>QC Metrics Summary Table</h2>
{summary_table}

<h2>Summary Charts</h2>
<div class="charts-grid">
{canvas_html}
</div>

<h2>QC Thresholds Reference</h2>
{threshold_table}

<h2>문제 해결 가이드 (Troubleshooting Guide)</h2>

<div class="note">
  <b>기본 진단 원칙:</b> Control Read 메트릭이 정상(≥500 reads, ≥50 kb, 일치율 ≥88%)이면
  문제의 원인은 <b>샘플 또는 라이브러리 준비</b>에 있습니다.
  Control이 비정상이면 <b>장비 또는 소모품</b>을 먼저 점검하십시오.
</div>

<!-- 진단 흐름도 -->
<div class="ts-section">
<h3>🔍 권장 진단 순서 (Recommended Troubleshooting Workflow)</h3>
<div class="flowchart">
  <div class="flow-step">
    <div class="flow-num">1</div>
    <div class="flow-body"><b>P1(%) 확인</b> — ZMW 로딩 수준 평가<br>
      최적 범위: 50–70%. 범위를 벗어나면 On-Plate Loading Concentration (OPLC)을 조정합니다.</div>
  </div>
  <div class="flow-step">
    <div class="flow-num">2</div>
    <div class="flow-body"><b>Control Read 메트릭 확인</b> — 장비·소모품 정상 여부 판단<br>
      Control reads ≥500, 평균 길이 ≥50 kb, 일치율 ≥88% 여부를 확인합니다.</div>
  </div>
  <div class="flow-step">
    <div class="flow-num">3</div>
    <div class="flow-body"><b>HiFi Yield 및 Read Length 확인</b> — 샘플/라이브러리 품질 평가<br>
      Control이 정상인데 Yield 또는 Read Length가 낮으면 gDNA 품질 또는 라이브러리 준비 문제입니다.</div>
  </div>
  <div class="flow-step">
    <div class="flow-num">4</div>
    <div class="flow-body"><b>Read Quality 및 Q30+ 확인</b> — 시퀀싱 화학 품질 평가<br>
      Median QV &lt;Q30이면 라이브러리 품질 또는 소모품 상태를 점검합니다.</div>
  </div>
  <div class="flow-step">
    <div class="flow-num">5</div>
    <div class="flow-body"><b>Missing Adapter 확인</b> — SMRTbell 라이브러리 완전성 평가<br>
      ≥10%이면 라이브러리 준비 단계의 어댑터 라이게이션 효율을 점검합니다.</div>
  </div>
  <div class="flow-step">
    <div class="flow-num">6</div>
    <div class="flow-body"><b>샘플 셋업 절차 재확인</b><br>
      위 항목들이 모두 정상이면 샘플 준비 및 SMRT Cell 셋업 절차를 재검토합니다.</div>
  </div>
</div>
</div>

<!-- P1 문제 -->
<div class="ts-section">
<h3>📊 P1 (ZMW Productivity) 이상</h3>
<div class="ts-grid">
  <div class="ts-card fail">
    <h4>❌ P1 &lt; 40% — Under-loading (과소 로딩)</h4>
    <ul>
      <li>빈 ZMW(P0)가 증가하여 데이터 생산량이 감소</li>
      <li>On-Plate Loading Concentration (OPLC)을 높여 재실험</li>
      <li>처음 배치 실험 시 Loading Titration으로 최적 OPLC 결정 권장</li>
      <li>목표 P1 범위: 55–65% (첫 실험 시), 최적화 후 50–70%</li>
    </ul>
    <div class="cause">원인: 샘플 농도 부족, 불충분한 어닐링, 라이브러리 희석 과다</div>
  </div>
  <div class="ts-card fail">
    <h4>❌ P1 &gt; 80% — Over-loading (과다 로딩)</h4>
    <ul>
      <li>Multi-loaded ZMW(P2) 증가 → 배경 잡음 상승, 데이터 품질 저하</li>
      <li>Mean Polymerase Read Length가 현저히 감소할 수 있음</li>
      <li>OPLC를 낮추어 재실험</li>
      <li>P2가 5% 초과 시 과다 로딩으로 판단</li>
    </ul>
    <div class="cause">원인: 샘플 농도 과다, OPLC 과다 설정</div>
  </div>
</div>
</div>

<!-- Control Read 문제 -->
<div class="ts-section">
<h3>🧪 Control Read 이상</h3>
<div class="ts-grid">
  <div class="ts-card fail">
    <h4>❌ Control Read Count &lt; 500</h4>
    <ul>
      <li>장비 또는 소모품 문제 가능성 높음 → 장비 점검 우선</li>
      <li>Revio Sequencing Control 시약 상태 확인 (유통기한, 보관 조건)</li>
      <li>SMRT Cell tray 및 Sequencing plate 상태 확인</li>
      <li>샘플 OPLC가 너무 높으면 샘플이 Control과 경쟁 → Control count 감소 가능</li>
    </ul>
    <div class="cause">원인: 소모품 불량, 장비 이상, 샘플 과다 로딩</div>
  </div>
  <div class="ts-card warn">
    <h4>⚠️ Control Read Length &lt; 50 kb (문서의 24-hr 기준선)</h4>
    <ul>
      <li><b>먼저 단위를 확인하십시오.</b> 이 값은 control 의 <b>polymerase read
          length</b>(한 ZMW 의 모든 통과 합) 여야 합니다. subread 길이 평균을
          쓰면 23~25 kb 로 나와서 항상 경고가 뜹니다 — 과거에 그 버그로 28샘플
          전부가 경고였습니다. 현재는 장비의 control.report.json 을 읽습니다.</li>
      <li>장비 및 소모품 성능 저하 의심</li>
      <li>Sequencing plate 개봉 후 사용 시간 확인 (Open Limit: 200시간)</li>
      <li>Polymerase kit 유통기한 및 보관 상태 점검</li>
      <li>Control Read Length ≥ 40 kb이면 샘플 문제에 집중하여 점검</li>
    </ul>
    <div class="cause">원인: 소모품 노후화, 부적절한 보관 (실온 노출, 빛 노출, 표면 오염)</div>
  </div>
</div>
</div>

<!-- HiFi Yield / Read Length 문제 -->
<div class="ts-section">
<h3>📉 HiFi Yield / Read Length 이상 (Control 정상인 경우)</h3>
<div class="ts-grid">
  <div class="ts-card warn">
    <h4>⚠️ HiFi Yield 낮음 (&lt;100 Gb, 15–20 kb WGS 행 기준)</h4>
    <ul>
      <li>Control이 정상이면 샘플 품질 또는 라이브러리 문제</li>
      <li>Insert size distribution 확인: Read Length Density plot 참조</li>
      <li>P1이 정상 범위(50–70%)인지 확인 — P1 낮으면 OPLC 조정이 우선</li>
      <li>15–20 kb library: 기대 Yield 100–120 Gb / SMRT Cell (문서 표는 24-hr 기준.
          30 hr 로 돌렸다면 이보다 유리한 조건이므로 100 Gb 미달은 더 분명한 신호다)</li>
      <li>Yield = P1 reads × Polymerase Read Length 임을 참고</li>
    </ul>
    <div class="cause">원인: 낮은 P1, 짧은 insert size, DNA 품질 불량, 로딩 농도 부적절</div>
  </div>
  <div class="ts-card warn">
    <h4>⚠️ HiFi Read Length 짧음 (&lt;15 kb)</h4>
    <ul>
      <li>라이브러리 insert size가 목표값보다 작음 → 라이브러리 준비 재검토</li>
      <li>gDNA 추출 과정에서 전단(shearing) 과다 여부 확인</li>
      <li>SMRTbell 라이브러리 준비 시 bead cleanup 조건 확인</li>
      <li>5 kb 미만의 minor peak는 잔존 short fragment로 정상 범위일 수 있음</li>
      <li>Bioanalyzer/Fragment Analyzer로 라이브러리 size distribution 재확인</li>
    </ul>
    <div class="cause">원인: gDNA 분해, 과도한 전단, 라이브러리 준비 조건 불량</div>
  </div>
</div>
</div>

<!-- Read Quality 문제 -->
<div class="ts-section">
<h3>🔬 Read Quality / Q30+ 이상</h3>
<div class="ts-grid">
  <div class="ts-card warn">
    <h4>⚠️ Median Read Quality &lt; Q30 또는 Q30+ Bases &lt; 90%</h4>
    <ul>
      <li>HiFi read는 QV ≥ 20이 최소 요건이며 일반적으로 QV 30–35 기대</li>
      <li>Mean passes가 낮으면(≤5) insert size에 비해 Polymerase Read Length가 짧은 것</li>
      <li>소모품(Polymerase kit, Sequencing plate) 상태 재점검</li>
      <li>Control concordance ≥ 88% 여부를 함께 확인</li>
    </ul>
    <div class="cause">원인: 낮은 subread pass 수, 소모품 품질 저하, DNA 손상</div>
  </div>
  <div class="ts-card warn">
    <h4>⚠️ Missing Adapter ≥ 10%</h4>
    <ul>
      <li>SMRTbell template에서 어댑터 1개만 검출된 ZMW 비율이 높음</li>
      <li>라이브러리 준비 시 어댑터 라이게이션 효율 확인</li>
      <li>SMRTbell prep kit 시약 상태 및 프로토콜 준수 여부 재검토</li>
      <li>&lt; 5%이면 정상, 5–10%이면 주의 모니터링</li>
    </ul>
    <div class="cause">원인: 라이게이션 효율 저하, 불완전한 SMRTbell 형성, 시약 불량</div>
  </div>
</div>
</div>

<!-- 활성 ZMW 감소 -->
<div class="ts-section">
<h3>📉 % Active ZMWs 급격한 감소 (Run Preview 확인 시)</h3>
<div class="ts-grid">
  <div class="ts-card info">
    <h4>ℹ️ Active ZMW 급감 (Run 중 모니터링)</h4>
    <ul>
      <li>최적 로딩 시 최대 % Active ZMWs ≈ 40% (P1 50–75% 해당)</li>
      <li>급격한 감소는 샘플 품질 문제 또는 소모품 불량 신호</li>
      <li>gDNA 과도한 손상 또는 Polymerase inhibitor 잔류 가능성</li>
      <li>bead cleanup 조건 재검토 (샘플 준비 단계부터 오염원 추적)</li>
      <li>Sequencing plate 보관 조건 확인 (실온 노출, 빛 노출 금지)</li>
    </ul>
    <div class="cause">원인: DNA 손상, Polymerase inhibitor, 소모품 부적절 보관</div>
  </div>
</div>
</div>

<!-- 참고 -->
<div class="note" style="margin-top:16px;">
  <b>추가 참고 문서:</b><br>
  • <i>Revio system run evaluation and troubleshooting guide</i> (PN 103-380-300) — 상세 트러블슈팅<br>
  • <i>SMRT Link Run Details report overview</i> (PN 103-604-700 Rev 01) — QC 메트릭 해석<br>
  • <i>Revio system specification sheet</i> (PN 102-326-552) — 최신 성능 사양
</div>

</div>
<footer>PacBio Raw QC Report &mdash; Generated by pacbio_rawdata_qc.py</footer>
{charts_js}
</body>
</html>"""


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="PacBio Revio raw data QC report generator"
    )
    parser.add_argument(
        "--run-dir", required=True, nargs="+",
        help="Path(s) to run directory. Multiple directories can be provided for multi-run projects."
    )
    parser.add_argument(
        "--sample-map", nargs="+", default=[],
        metavar="WELL=SAMPLE",
        help="Map wells to sample names, e.g. 1_A01=TG-6102 1_B01=TG-6283"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output HTML path (default: QC_Report_rawdata.html in first run-dir)"
    )
    parser.add_argument(
        "--movie-time", type=int, default=30,
        help="이 run 의 실제 movie collection time (시간). 표시용이며 임계값을 "
             "바꾸지 않는다. 우리 장비의 6개 run 은 모두 30 hr 이었다 "
             "(MovieLength=1800.0 분). default: 30"
    )
    parser.add_argument(
        "--genome-size", type=float, default=2.7,
        help="Genome size in Gb for coverage estimation (default: 2.7 for mouse GRCm39)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run_dirs = [d.rstrip("/") for d in args.run_dir]

    # sample mapping
    sample_map = {}
    for item in args.sample_map:
        if "=" in item:
            well, name = item.split("=", 1)
            sample_map[well] = name

    all_metrics = []
    for run_dir in run_dirs:
        wells = sorted([
            d for d in os.listdir(run_dir)
            if os.path.isdir(os.path.join(run_dir, d)) and re.match(r"\d_[A-D]\d+", d)
        ])
        if not wells:
            print(f"[WARN] No well directories found in {run_dir}", file=sys.stderr)
            continue
        print(f"Run: {os.path.basename(run_dir)}  →  Wells: {wells}")
        for well in wells:
            sample = sample_map.get(well, well)
            print(f"  Processing {well} ({sample}) ...")
            m = collect_well_metrics(run_dir, well, sample, args.genome_size)
            if m:
                all_metrics.append(m)

    if not all_metrics:
        print("ERROR: No metrics collected.", file=sys.stderr)
        sys.exit(1)

    output = args.output or os.path.join(run_dirs[0], "QC_Report_rawdata.html")
    html = build_html(all_metrics, run_dirs, args.movie_time, args.genome_size)
    with open(output, "w") as f:
        f.write(html)

    print(f"\nReport saved: {output}")
    print("\nOverall results:")
    for m in all_metrics:
        print(f"  [{m['run']}] {m['well']} ({m['sample']}) [{m['barcode']}]: {m['overall']}")


if __name__ == "__main__":
    main()
