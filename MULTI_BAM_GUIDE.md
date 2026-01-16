# HiFi 데이터 파일(BAM)이 여러 개일 때 파이프라인 구성 가이드

HiFi-human-WGS-WDL 워크플로우에서 여러 개의 BAM 파일을 처리하는 방법은 **사용 사례**에 따라 달라집니다.

## 시나리오별 구성 방법

### 1. 한 개인(Sample)의 여러 BAM 파일 → `singleton` 워크플로우

**상황**: 동일한 개인의 HiFi 데이터가 여러 BAM 파일로 나뉘어 있는 경우 (예: 여러 SMRT Cell에서 시퀀싱)

**해결책**: `hifi_reads` 배열에 모든 BAM 파일 경로를 나열

#### 입력 파일 예시 (`singleton.inputs.json`):
```json
{
  "humanwgs_singleton.sample_id": "SAMPLE001",
  "humanwgs_singleton.sex": "MALE",
  "humanwgs_singleton.hifi_reads": [
    "/path/to/data/sample001_cell1.bam",
    "/path/to/data/sample001_cell2.bam",
    "/path/to/data/sample001_cell3.bam"
  ],
  "humanwgs_singleton.ref_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.ref_map.v3p1p0.hpc.tsv",
  "humanwgs_singleton.tertiary_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.tertiary_map.v3p1p0.hpc.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": true
}
```

#### 실행 명령:
```bash
conda activate hifi-human-wgs
miniwdl run workflows/singleton.wdl --input singleton.inputs.json
```

**동작 방식**: 워크플로우가 각 BAM 파일을 개별적으로 정렬(align)한 후 자동으로 병합(merge)합니다.

---

### 2. 여러 개인(Samples)의 BAM 파일 → `family` 워크플로우

**상황**: 여러 사람의 HiFi 데이터를 함께 분석하려는 경우 (가족, 코호트 등)

**해결책**: `family` 워크플로우 사용하여 joint calling 수행

#### 입력 파일 예시 (`family.inputs.json`):
```json
{
  "humanwgs_family.family": {
    "family_id": "FAM001",
    "samples": [
      {
        "sample_id": "FATHER",
        "hifi_reads": [
          "/path/to/data/father_cell1.bam",
          "/path/to/data/father_cell2.bam"
        ],
        "affected": false,
        "sex": "MALE",
        "father_id": null,
        "mother_id": null
      },
      {
        "sample_id": "MOTHER",
        "hifi_reads": [
          "/path/to/data/mother_cell1.bam"
        ],
        "affected": false,
        "sex": "FEMALE",
        "father_id": null,
        "mother_id": null
      },
      {
        "sample_id": "CHILD",
        "hifi_reads": [
          "/path/to/data/child_cell1.bam",
          "/path/to/data/child_cell2.bam"
        ],
        "affected": true,
        "sex": "MALE",
        "father_id": "FATHER",
        "mother_id": "MOTHER"
      }
    ]
  },
  "humanwgs_family.phenotypes": "HP:0001250,HP:0001263",
  "humanwgs_family.ref_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.ref_map.v3p1p0.hpc.tsv",
  "humanwgs_family.tertiary_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.tertiary_map.v3p1p0.hpc.tsv",
  "humanwgs_family.backend": "HPC",
  "humanwgs_family.preemptible": true
}
```

#### 실행 명령:
```bash
conda activate hifi-human-wgs
miniwdl run workflows/family.wdl --input family.inputs.json
```

**장점**:
- Joint variant calling으로 가족 내 변이 정확도 향상
- 멘델 유전 패턴 분석 가능
- De novo 변이 검출 가능

---

### 3. 여러 독립적인 개인의 BAM 파일 → 여러 번 실행 또는 배치 처리

**상황**: 관련 없는 여러 개인의 샘플을 각각 독립적으로 분석

**옵션 A: 각 샘플에 대해 워크플로우를 개별 실행**

```bash
# Sample 1 실행
miniwdl run workflows/singleton.wdl --input sample1.inputs.json

# Sample 2 실행
miniwdl run workflows/singleton.wdl --input sample2.inputs.json

# Sample 3 실행
miniwdl run workflows/singleton.wdl --input sample3.inputs.json
```

**옵션 B: 배치 스크립트 작성**

프로젝트에 포함된 `batch_run.sh` 스크립트를 사용하여 여러 샘플을 효율적으로 처리할 수 있습니다.

**1. 스크립트 준비**
```bash
# 실행 권한 부여
chmod +x batch_run.sh

# 입력 파일 디렉토리 구조 생성
mkdir -p inputs
mkdir -p outputs
mkdir -p logs
```

**2. 입력 파일 준비**

각 샘플의 입력 JSON 파일을 `inputs/` 디렉토리에 생성:
```bash
# inputs/sample1.inputs.json
# inputs/sample2.inputs.json
# inputs/sample3.inputs.json
# 등...
```

**3. 샘플 목록 설정**

`batch_run.sh` 파일을 편집하여 SAMPLES 배열을 수정:
```bash
nano batch_run.sh

# 다음 부분을 실제 샘플 이름으로 수정
SAMPLES=("sample1" "sample2" "sample3" "sample4")
```

**4. 스크립트 실행**

```bash
# 병렬 실행 (권장 - 여러 샘플 동시 처리)
./batch_run.sh parallel

# 순차 실행 (한 번에 하나씩 처리)
./batch_run.sh sequential
```

**5. 진행 상황 확인**

```bash
# 로그 파일 확인
tail -f logs/sample1.log

# 모든 로그 파일 확인
ls -lh logs/

# 특정 샘플의 출력 디렉토리 확인
ls -lh outputs/sample1/
```

**배치 스크립트의 주요 기능:**
- ✅ 병렬/순차 실행 모드 선택 가능
- ✅ 자동 디렉토리 생성 (outputs, logs)
- ✅ 입력 파일 존재 여부 자동 검증
- ✅ 각 샘플별 독립적인 로그 파일 생성
- ✅ 실행 시간 측정 및 요약 정보 제공
- ✅ Conda 환경 활성화 확인

**간단한 배치 스크립트 예시 (수동 작성):**
```bash
#!/bin/bash
# simple_batch.sh

SAMPLES=("sample1" "sample2" "sample3" "sample4")

for sample in "${SAMPLES[@]}"; do
    echo "Processing $sample..."
    miniwdl run workflows/singleton.wdl \
        --input inputs/${sample}.inputs.json \
        --dir outputs/${sample} \
        > logs/${sample}.log 2>&1 &
done

wait
echo "All samples processed!"
```

**옵션 C: 병렬 처리 (GNU Parallel 사용)**

```bash
# inputs 디렉토리의 모든 입력 파일 처리
ls inputs/*.inputs.json | parallel -j 4 \
    'miniwdl run workflows/singleton.wdl --input {} --dir outputs/{/.}'
```

---

### 4. Fail Reads가 있는 경우

일부 샘플에서 품질이 낮은 reads(fail reads)를 별도로 처리하고 싶을 때:

```json
{
  "humanwgs_singleton.sample_id": "SAMPLE001",
  "humanwgs_singleton.hifi_reads": [
    "/path/to/data/hifi_pass_reads.bam"
  ],
  "humanwgs_singleton.fail_reads": [
    "/path/to/data/hifi_fail_reads_cell1.bam",
    "/path/to/data/hifi_fail_reads_cell2.bam"
  ],
  "humanwgs_singleton.sex": "FEMALE",
  "humanwgs_singleton.ref_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.ref_map.v3p1p0.hpc.tsv",
  "humanwgs_singleton.tertiary_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.tertiary_map.v3p1p0.hpc.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": true
}
```

**참고**: Fail reads는 TRGT (Tandem Repeat Genotyping) 영역에만 매핑되어 해당 영역의 커버리지를 향상시킵니다.

---

## 실전 예제

### 예제 1: 단일 샘플, 3개의 SMRT Cell BAM 파일

```bash
# 1. 입력 파일 준비
cat > my_sample.inputs.json << 'EOF'
{
  "humanwgs_singleton.sample_id": "Patient_A",
  "humanwgs_singleton.sex": "MALE",
  "humanwgs_singleton.hifi_reads": [
    "/data/hifi/patient_a/m64012_201201_012345.hifi_reads.bam",
    "/data/hifi/patient_a/m64012_201202_123456.hifi_reads.bam",
    "/data/hifi/patient_a/m64012_201203_234567.hifi_reads.bam"
  ],
  "humanwgs_singleton.ref_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.ref_map.v3p1p0.hpc.tsv",
  "humanwgs_singleton.tertiary_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.tertiary_map.v3p1p0.hpc.tsv",
  "humanwgs_singleton.backend": "HPC",
  "humanwgs_singleton.preemptible": true
}
EOF

# 2. 워크플로우 실행
conda activate hifi-human-wgs
miniwdl run workflows/singleton.wdl --input my_sample.inputs.json --verbose
```

### 예제 2: Trio 분석 (부모-자녀)

```bash
cat > trio_family.inputs.json << 'EOF'
{
  "humanwgs_family.family": {
    "family_id": "TRIO_001",
    "samples": [
      {
        "sample_id": "PROBAND",
        "hifi_reads": ["/data/proband.bam"],
        "affected": true,
        "sex": "MALE",
        "father_id": "FATHER",
        "mother_id": "MOTHER"
      },
      {
        "sample_id": "FATHER",
        "hifi_reads": ["/data/father.bam"],
        "affected": false,
        "sex": "MALE"
      },
      {
        "sample_id": "MOTHER",
        "hifi_reads": ["/data/mother.bam"],
        "affected": false,
        "sex": "FEMALE"
      }
    ]
  },
  "humanwgs_family.phenotypes": "HP:0001263,HP:0001250",
  "humanwgs_family.ref_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.ref_map.v3p1p0.hpc.tsv",
  "humanwgs_family.tertiary_map_file": "/home/ygkim/hifi-wdl-resources/dataset/GRCh38.tertiary_map.v3p1p0.hpc.tsv",
  "humanwgs_family.backend": "HPC",
  "humanwgs_family.preemptible": true
}
EOF

miniwdl run workflows/family.wdl --input trio_family.inputs.json
```

---

## 파이프라인 선택 가이드

| 상황 | 사용할 워크플로우 | 주요 특징 |
|------|------------------|----------|
| 한 사람의 여러 BAM 파일 | `singleton.wdl` | `hifi_reads` 배열에 모든 파일 나열 |
| 가족/Trio 분석 | `family.wdl` | Joint calling, de novo 변이 검출 |
| 여러 독립 샘플 | `singleton.wdl` (각각) | 배치 처리 또는 순차 실행 |
| 대규모 코호트 (10+ 샘플) | `family.wdl` 또는 개별 `singleton` + joint analysis | 용도에 따라 선택 |

---

## 추가 팁

### 1. BAM 파일 병합을 미리 하는 것이 좋을까?
**권장하지 않음**. 워크플로우가 자동으로 병합하므로 원본 BAM을 그대로 사용하는 것이 좋습니다.

### 2. 리소스 요구사항
- 각 샘플당 최대 64 CPU cores, 256 GB RAM 필요
- 디스크 공간: 샘플당 약 500GB-1TB (중간 파일 포함)

### 3. 입력 파일 경로
- 절대 경로 사용 권장: `/full/path/to/file.bam`
- 상대 경로는 워크플로우 실행 위치에 따라 달라질 수 있음

### 4. 샘플 ID 명명 규칙
- 영문자, 숫자, 마침표(.), 하이픈(-), 언더스코어(_)만 사용
- 공백이나 특수문자는 피할 것

---

## 문제 해결

### Q: "File not found" 오류가 발생합니다
**A**: 
- BAM 파일 경로가 올바른지 확인
- 절대 경로 사용
- 파일 접근 권한 확인

### Q: 메모리 부족 오류
**A**: 
- HPC 설정에서 메모리 할당량 증가
- `miniwdl.cfg`에서 리소스 조정

### Q: 여러 샘플을 동시에 실행하고 싶습니다
**A**: 
- 배치 스크립트 또는 GNU Parallel 사용
- SLURM의 경우 job array 기능 활용

---

## 참고 문서

- [Singleton 워크플로우 상세 문서](docs/singleton.md)
- [Family 워크플로우 상세 문서](docs/family.md)
- [HPC 백엔드 설정](docs/backend-hpc.md)
- [입력 파일 참조 맵](docs/ref_map.md)
