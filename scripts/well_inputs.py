#!/usr/bin/env python3
"""웰 metadata 에서 humanwgs_singleton 입력 json 을 만든다.

왜 이게 있나: 샘플 이름·바코드·movie 이름을 **손으로 적다가 틀린 적이 두 번
있다** (P8, 낡은 입력 경로). 세 값 모두 웰 디렉터리가 스스로 말한다 —

    <well>/metadata/*.metadata.xml   의 <WellSample Name="...">  ← 샘플 이름
    <well>/hifi_reads/*.hifi_reads.bc####.bam                    ← 바코드·경로

그래서 읽어서 쓴다. BioSample Name 은 쓰지 않는다 — 그건 별도 라벨이고
(예: 2026-11-kyg-01) 분석에서 쓰는 이름이 아니다.

사용:
    python3 scripts/well_inputs.py <run_dir> <well> [<well> ...]
    python3 scripts/well_inputs.py <run_dir> <well> --write   # 파일로 저장

전송이 끝나지 않은 웰은 건너뛰고 이유를 말한다. 전송 중 임시 파일
(`.<name>.bam.<suffix>`) 을 완성된 BAM 으로 착각하지 않는다.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import xml.etree.ElementTree as ET

# 이 저장소의 고정 설정. 바꾸려면 여기서 바꾸고 커밋한다.
REF_MAP = "/data_4tb/hifi-human-wgs-wdl-custom/GRCm39.ref_map.tsv"
DEFAULTS = {
    "backend": "HPC",
    "preemptible": False,
    "gpu": True,
    "run_pgx": False,
    "run_assembly": True,
    "busco_lineage": "mammalia_odb10",
}
WF = "humanwgs_singleton"


def well_sample_name(well_dir: str) -> str | None:
    """`WellSample Name` 을 읽는다. 없으면 None."""
    xmls = sorted(glob.glob(os.path.join(well_dir, "metadata", "*.metadata.xml")))
    if not xmls:
        return None
    # full_preview 쪽이 더 늦게 쓰이므로 있으면 그걸 쓴다.
    xmls.sort(key=lambda p: ("full_preview" in p, p))
    for path in reversed(xmls):
        try:
            root = ET.parse(path).getroot()
        except ET.ParseError:
            continue
        for el in root.iter():
            # 네임스페이스가 붙으므로 태그 끝만 본다.
            if el.tag.rsplit("}", 1)[-1] == "WellSample":
                name = el.attrib.get("Name")
                if name:
                    return name
    return None


def hifi_bam(well_dir: str) -> tuple[str | None, str]:
    """완성된 hifi_reads BAM 경로와 상태 메시지를 돌려준다."""
    d = os.path.join(well_dir, "hifi_reads")
    if not os.path.isdir(d):
        return None, "hifi_reads 디렉터리가 없다 (전송 시작 전)"
    bams = [
        p
        for p in glob.glob(os.path.join(d, "*.hifi_reads.bc*.bam"))
        if not os.path.basename(p).startswith(".")
    ]
    if not bams:
        # 전송 중 임시 파일이 있으면 그 사실을 말해 준다.
        partial = glob.glob(os.path.join(d, ".*.bam.*"))
        if partial:
            size = sum(os.path.getsize(p) for p in partial) / 1e9
            return None, f"전송 중이다 ({size:.1f} GB, 임시 파일 {len(partial)}개)"
        return None, "완성된 BAM 이 없다"
    if len(bams) > 1:
        return None, f"BAM 이 {len(bams)}개다 — 어느 것인지 사람이 정해야 한다"
    return bams[0], "ok"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("wells", nargs="+")
    ap.add_argument(
        "--write",
        action="store_true",
        help="<sample>.inputs.json 으로 저장한다 (기존 파일은 덮지 않는다)",
    )
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    rc = 0
    for well in args.wells:
        well_dir = os.path.join(args.run_dir, well)
        if not os.path.isdir(well_dir):
            print(f"[{well}] 건너뜀: {well_dir} 가 없다", file=sys.stderr)
            rc = 1
            continue

        name = well_sample_name(well_dir)
        bam, status = hifi_bam(well_dir)

        if name is None:
            print(f"[{well}] 건너뜀: metadata.xml 이 아직 없어 샘플 이름을 모른다",
                  file=sys.stderr)
            rc = 1
            continue
        if bam is None:
            print(f"[{well}] 건너뜀: {status}  (샘플 이름은 '{name}')", file=sys.stderr)
            rc = 1
            continue

        bc = re.search(r"\.(bc\d+)\.bam$", os.path.basename(bam))
        inputs = {
            f"{WF}.sample_id": name,
            f"{WF}.hifi_reads": [bam],
            f"{WF}.ref_map_file": REF_MAP,
            **{f"{WF}.{k}": v for k, v in DEFAULTS.items()},
        }
        body = json.dumps(inputs, indent=2, ensure_ascii=False) + "\n"
        print(f"[{well}] sample_id='{name}' barcode={bc.group(1) if bc else '?'}",
              file=sys.stderr)

        if args.write:
            out = os.path.join(args.out_dir, f"{name}.inputs.json")
            if os.path.exists(out):
                print(f"[{well}] {out} 가 이미 있다 — 덮지 않는다", file=sys.stderr)
                rc = 1
                continue
            with open(out, "w") as fh:
                fh.write(body)
            print(f"[{well}] -> {out}  (git add -f 로 추적에 넣는다)", file=sys.stderr)
        else:
            print(body)
    return rc


if __name__ == "__main__":
    sys.exit(main())
