# Base Modification Discovery — HPC 전환 계획

## 배경

파일럿 실행(그룹당 1샘플, 총 4샘플)을 로컬 서버에서 실행한 결과, jasmine 단계만으로 샘플당 약 40~55시간 소요 확인.
전체 12샘플을 로컬 순차 실행할 경우 현실적인 완료가 불가능하므로 HPC 클러스터 전환 필요.

---

## 현재 로컬 실행 환경

| 항목 | 내용 |
|---|---|
| WDL runner | miniwdl (singularity 백엔드) |
| 설정 파일 | `/data_4tb/hifi-human-wgs-wdl-custom/config/miniwdl.local.cfg` |
| `task_concurrency` | **1** (1개 task만 순차 실행) |
| 입력 BAM 위치 | `/mnt/JJ_dis_8tb/r84285_20260326_071034/` (NTFS USB 외장 HDD, ntfs-3g 마운트) |
| 출력 위치 | `/mnt/JJ_dis_8tb/base_mod_discovery_pilot_output/` |
| Call cache | `/data_4tb/hifi-human-wgs-wdl-custom/miniwdl_call_cache/` |

### 로컬 실행 시 소요 시간 추정 (12샘플 전체)

| Task | 건수 | 샘플당 시간 | 합계 (순차) |
|---|---|---|---|
| pbmm2_align_kinetics | 12 | ~6h | ~72h |
| jasmine_modification_calling | 12 | ~48h | **~576h (24일)** |
| ipd_summary | 4 비교 | ~수시간 | TBD |
| jasmine_pileup (modkit) | 12 | ~2h | ~24h |
| 이후 downstream | - | 빠름 | <수시간 |

→ **jasmine 단독으로 24일, 현실적으로 불가**

---

## HPC 전환 목표

- WDL `scatter` 병렬 실행 활용: 동일 step 내 모든 샘플 동시 처리
- jasmine 12샘플 동시 제출 시 wall time ≈ 최대 샘플 1개 분량 (~55h)
- 전체 파이프라인 완료 목표: **약 3~5일 (클러스터 큐 대기 포함)**

---

## 전환 체크리스트

### 1. HPC 환경 확인
- [ ] 클러스터 스케줄러 종류 확인 (SLURM / PBS / SGE)
- [ ] 노드당 CPU/메모리 사양 확인 (jasmine: GPU 권장, 최소 32코어/64GB)
- [ ] Singularity / Apptainer 설치 여부 확인
- [ ] 인터넷 접근 가능 여부 (컨테이너 이미지 pull 가능한지)

### 2. 데이터 이전
- [ ] 입력 BAM 12개를 클러스터 스토리지로 복사
  - 현재: `/mnt/JJ_dis_8tb/r84285_20260326_071034/1_{A,B,C,D}0{1..3}/hifi_reads/*.bam`
  - 대상: 클러스터 공유 스토리지 (예: `/scratch/` 또는 프로젝트 디렉토리)
- [ ] 레퍼런스 파일 복사 또는 마운트 확인
  - GRCm39 FASTA: `/data_4tb/hifi-human-wgs-wdl-custom/hifi-wdl-resources/GRCm39/`
  - Gencode vM36 GTF, CpG island TSV 포함
- [ ] Singularity 이미지 캐시 이전 (인터넷 불가 클러스터인 경우)
  - 스크립트: `./scripts/populate_miniwdl_singularity_cache.sh`

### 3. miniwdl HPC 설정 (SLURM 기준)

```bash
# miniwdl-slurm 설치
pip install miniwdl-slurm
```

`~/.config/miniwdl.cfg` (클러스터용):
```ini
[scheduler]
task_concurrency = 200       # 클러스터는 높게 설정

[singularity]
# Singularity 이미지 캐시 경로 (공유 스토리지)
image_cache = /scratch/ygkim/miniwdl_singularity_cache

[file_io]
allow_any_input = true        # map 파일 접근 허용
output_hardlinks = false

[SLURM]                       # miniwdl-slurm 섹션
# 기본 파티션/큐 설정 (클러스터별 조정 필요)
# 참고: https://github.com/miniwdl-ext/miniwdl-slurm
```

> 참고 문서: `docs/backend-hpc.md` — miniwdl-slurm 및 Cromwell 설정 가이드 포함

### 4. inputs.json 수정

```json
"base_modification_discovery.control_bams": [
  "/scratch/ygkim/r84285/1_A01/hifi_reads/...bam",  ← 클러스터 경로로 변경
  "/scratch/ygkim/r84285/1_A02/hifi_reads/...bam",
  "/scratch/ygkim/r84285/1_A03/hifi_reads/...bam"   ← 3샘플로 확장
],
...
"base_modification_discovery.default_runtime_attributes": {
  "backend": "HPC"   ← 이미 설정됨
}
```

### 5. 파일럿 결과 재사용 (call cache)

파일럿에서 완료된 task는 call cache 활용 가능:
- `pbmm2_align_kinetics` × 4샘플 ✓
- `jasmine_modification_calling` × 4샘플 ✓ (완료 예정)

단, call cache는 로컬 경로 기반이므로 **클러스터에서는 새로 실행 필요**.
파일럿 분석 결과(BAM, BED 등)는 수동으로 output 디렉토리에서 가져와 재사용 검토 가능.

---

## 대안: 로컬 부분 완료 후 HPC 전환 전략

파일럿(4샘플)은 현재 로컬에서 완료까지 진행 → 결과 검증 후 전체 12샘플은 HPC에서 실행.

| 단계 | 환경 | 목적 |
|---|---|---|
| 파일럿 4샘플 전체 파이프라인 | 로컬 (현재 진행 중) | 파이프라인 검증, 파라미터 최적화 |
| 전체 12샘플 | **HPC** | 실제 분석 |

---

## 참고

- 현재 파일럿 로그: `/data_4tb/hifi-human-wgs-wdl-custom/output/base_mod_discovery_pilot.nohup.log`
- 파이프라인 설계 문서: `docs/base_modification_discovery_design.md`
- HPC 백엔드 일반 설정: `docs/backend-hpc.md`
- inputs 템플릿: `base_mod_discovery.inputs.template.json`
