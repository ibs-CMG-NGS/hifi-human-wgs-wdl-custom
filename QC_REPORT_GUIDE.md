# HiFi WGS Pipeline - QC Report 가이드

## 📊 개요

HiFi-human-WGS-WDL 파이프라인은 배치 처리 완료 후 자동으로 **종합 QC 리포트(HTML)**를 생성하는 기능을 제공합니다. 이 리포트는 여러 샘플의 분석 결과를 시각적으로 요약하여 품질 관리와 결과 검토를 용이하게 합니다.

---

## 🎯 주요 기능

### 자동 생성되는 QC Metrics

1. **Read Statistics**
   - Total reads (총 리드 수)
   - Mean read length (평균 리드 길이)
   - Mean read quality (평균 리드 품질)
   - Read length distribution

2. **Alignment Statistics**
   - Mapping rate (매핑률)
   - Aligned bases
   - Primary/secondary alignments

3. **Coverage Statistics**
   - Mean coverage depth (평균 커버리지)
   - Coverage uniformity
   - Chromosome별 coverage

4. **Variant Calling Results**
   - Small variants (SNPs, Indels)
     - Total count
     - Het/Hom ratio
     - Ti/Tv ratio
   - Structural variants (SVs)
     - Deletions, Insertions, Duplications
     - Inversions, Translocations

5. **Phasing Statistics**
   - Phase block N50
   - Phased variants
   - Switch error rate

6. **PharmCAT Results**
   - 약물유전체 분석 결과
   - Phenotype 예측
   - 약물 반응성 정보

7. **File Information**
   - 주요 출력 파일 크기
   - 파일 생성 상태 체크
   - 디스크 사용량

---

## 🚀 사용 방법

### 1. 자동 생성 (권장)

배치 스크립트를 실행하면 완료 시 자동으로 QC 리포트가 생성됩니다:

```bash
# 병렬 실행 (자동으로 QC 리포트 생성)
./batch_run.sh

# 또는 최적화 버전
./batch_run_optimized.sh parallel
```

**생성 위치:**
```
/data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_YYYYMMDD_HHMMSS.html
```

### 2. 수동 생성

배치 작업 완료 후 별도로 리포트를 생성할 수도 있습니다:

```bash
# 모든 샘플 포함
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/hifi-human-wgs-wdl-custom/batch_results \
  --output QC_Report.html

# 특정 샘플만 포함
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/hifi-human-wgs-wdl-custom/batch_results \
  --output Custom_Report.html \
  --samples KTY9537 KTY9538 KTY9539

# 커스텀 출력 경로
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/hifi-human-wgs-wdl-custom/batch_results \
  --output /path/to/reports/MyReport_$(date +%Y%m%d).html
```

---

## 📖 명령줄 옵션

### `generate_qc_report.py` 파라미터

| 파라미터 | 설명 | 기본값 | 필수 |
|---------|------|-------|-----|
| `--batch-results` | 배치 결과 디렉토리 경로 | `/data_4tb/hifi-human-wgs-wdl-custom/batch_results` | 아니오 |
| `--output` | 출력 HTML 파일 경로 | `QC_Report.html` | 아니오 |
| `--samples` | 포함할 샘플 ID 목록 (공백으로 구분) | 모든 샘플 | 아니오 |

### 예제

```bash
# 1. 기본 사용
python3 scripts/generate_qc_report.py

# 2. 커스텀 출력 경로
python3 scripts/generate_qc_report.py \
  --output /mnt/reports/weekly_report.html

# 3. 특정 샘플 2개만
python3 scripts/generate_qc_report.py \
  --samples sample1 sample2

# 4. 다른 배치 디렉토리 지정
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/project_X/results \
  --output /data_4tb/project_X/reports/QC_Report.html
```

---

## 🔍 리포트 확인 방법

### Linux/WSL 환경

```bash
# 브라우저로 열기
firefox /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_*.html

# 또는 Google Chrome
google-chrome /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_*.html

# 최신 리포트 자동 열기
firefox $(ls -t /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_*.html | head -1)
```

### Windows에서 WSL 파일 접근

WSL Ubuntu 파일 시스템은 Windows에서 다음 경로로 접근 가능합니다:

```
\\wsl.localhost\Ubuntu\data_4tb\hifi-human-wgs-wdl-custom\batch_results\QC_Report_*.html
```

**방법 1: 파일 탐색기**
1. Windows 파일 탐색기 열기
2. 주소창에 `\\wsl.localhost\Ubuntu\data_4tb\hifi-human-wgs-wdl-custom\batch_results` 입력
3. QC_Report_*.html 파일을 더블클릭

**방법 2: PowerShell에서**
```powershell
# 최신 리포트 열기
Start-Process "\\wsl.localhost\Ubuntu\data_4tb\hifi-human-wgs-wdl-custom\batch_results\QC_Report_*.html"
```

---

## 📁 필요한 입력 파일

QC 리포트 생성을 위해 다음 파일들이 배치 결과 디렉토리에 있어야 합니다:

```
batch_results/
├── sample1/
│   ├── _LAST/                          # 최신 워크플로우 심볼릭 링크
│   │   ├── outputs.json                # ✅ 필수: 워크플로우 출력 정보
│   │   └── out/
│   │       ├── bam_statistics/
│   │       │   └── *_stats.txt         # ✅ 필수: Read/alignment 통계
│   │       └── mosdepth_summary/
│   │           └── *.mosdepth.summary.txt  # ✅ 필수: Coverage 통계
│   └── 20260130_123456_humanwgs_singleton/
│       └── (실제 워크플로우 디렉토리)
├── sample2/
│   └── ...
└── logs/
    ├── sample1.log
    └── sample2.log
```

### 파일 설명

1. **outputs.json**
   - 워크플로우의 모든 출력 파일 경로
   - Variant calling 결과 파일 정보
   - PharmCAT, TRGT 등 tertiary 분석 결과

2. **\*_stats.txt** (BAM statistics)
   - `samtools stats` 결과
   - Read length, quality, mapping 정보

3. **\*.mosdepth.summary.txt**
   - `mosdepth` 결과
   - Chromosome별 coverage 통계

---

## 📊 리포트 내용 상세

### 1. Summary Table (요약 테이블)

모든 샘플의 주요 QC metrics를 한눈에 확인:

| Sample ID | Mean Coverage | Total Reads | Mean Read Length | Mapping Rate | SNPs | Indels | SVs |
|-----------|---------------|-------------|------------------|--------------|------|--------|-----|
| KTY9537   | 34.2x        | 12,345,678  | 15,234 bp       | 99.2%        | 3.5M | 450K   | 12K |
| KTY9538   | 32.8x        | 11,987,654  | 14,987 bp       | 99.1%        | 3.4M | 445K   | 11K |

### 2. Individual Sample Reports

각 샘플별 상세 정보:

- **Read Quality Metrics**
  - Read length histogram
  - Base quality distribution
  - GC content

- **Alignment Metrics**
  - Mapping quality scores
  - Insert size distribution
  - Properly paired reads %

- **Coverage Analysis**
  - Genome-wide coverage plot
  - Coverage distribution
  - Low coverage regions

- **Variant Statistics**
  - Variant type breakdown (SNP/Indel/SV)
  - Quality score distribution
  - Depth distribution

### 3. Comparative Analysis

여러 샘플 비교:

- Coverage uniformity across samples
- Variant calling consistency
- Batch effect 탐지

### 4. QC Pass/Fail Criteria

자동 품질 평가:

| Metric | Threshold | Status |
|--------|-----------|--------|
| Mean Coverage | ≥ 30x | ✅ PASS |
| Mapping Rate | ≥ 95% | ✅ PASS |
| Mean Read Quality | ≥ Q20 | ✅ PASS |
| Mean Read Length | ≥ 10kb | ✅ PASS |

---

## 🛠️ 고급 사용법

### 1. 리포트 커스터마이징

스크립트를 수정하여 추가 metrics 포함:

```python
# scripts/generate_qc_report.py 수정

def parse_custom_metrics(sample_dir):
    """사용자 정의 metrics 추가"""
    # 예: VCF 파일에서 추가 통계 추출
    vcf_file = find_vcf(sample_dir)
    custom_stats = analyze_vcf(vcf_file)
    return custom_stats
```

### 2. 여러 배치 결과 통합

```bash
# 여러 배치의 결과를 하나의 리포트로
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/batch1 \
  --output combined_report.html

# 수동으로 샘플 추가 (심볼릭 링크 활용)
mkdir -p /data_4tb/combined_results
ln -s /data_4tb/batch1/sample1 /data_4tb/combined_results/
ln -s /data_4tb/batch2/sample2 /data_4tb/combined_results/

python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/combined_results \
  --output combined_report.html
```

### 3. 자동화된 리포트 이메일 전송

```bash
#!/bin/bash
# 리포트 생성 및 이메일 전송 자동화

REPORT="/data_4tb/batch_results/QC_Report_$(date +%Y%m%d).html"

python3 scripts/generate_qc_report.py --output "${REPORT}"

# 이메일 전송 (mailx 필요)
echo "QC Report 생성 완료" | mail \
  -s "HiFi WGS QC Report - $(date +%Y-%m-%d)" \
  -a "${REPORT}" \
  your.email@example.com
```

### 4. 리포트를 웹 서버에 게시

```bash
# 웹 서버 디렉토리로 복사
REPORT="/data_4tb/batch_results/QC_Report_$(date +%Y%m%d).html"
WEB_DIR="/var/www/html/qc_reports"

python3 scripts/generate_qc_report.py --output "${REPORT}"
sudo cp "${REPORT}" "${WEB_DIR}/"
sudo chmod 644 "${WEB_DIR}/$(basename ${REPORT})"

echo "Report available at: http://your-server/qc_reports/$(basename ${REPORT})"
```

---

## 🐛 트러블슈팅

### 문제 1: 리포트가 생성되지 않음

```bash
# 스크립트 존재 확인
ls -l scripts/generate_qc_report.py

# 수동 실행으로 에러 확인
python3 scripts/generate_qc_report.py \
  --batch-results /data_4tb/hifi-human-wgs-wdl-custom/batch_results \
  --output test_report.html
```

**가능한 원인:**
- Python 환경 문제
- 필수 파일(outputs.json, stats.txt) 누락
- 디렉토리 권한 문제

### 문제 2: 일부 샘플이 리포트에 포함되지 않음

```bash
# 샘플 디렉토리 구조 확인
ls -la /data_4tb/hifi-human-wgs-wdl-custom/batch_results/sample_name/

# _LAST 심볼릭 링크 확인
ls -la /data_4tb/hifi-human-wgs-wdl-custom/batch_results/sample_name/_LAST

# outputs.json 존재 확인
cat /data_4tb/hifi-human-wgs-wdl-custom/batch_results/sample_name/_LAST/outputs.json
```

**해결책:**
```bash
# _LAST 심볼릭 링크 수동 생성
cd /data_4tb/hifi-human-wgs-wdl-custom/batch_results/sample_name
ln -sf 20260130_123456_humanwgs_singleton _LAST
```

### 문제 3: Coverage 정보가 0으로 표시됨

```bash
# mosdepth 결과 파일 확인
find /data_4tb/hifi-human-wgs-wdl-custom/batch_results/sample_name \
  -name "*.mosdepth.summary.txt"

# 파일 내용 확인
cat /path/to/sample.mosdepth.summary.txt
```

**원인:**
- mosdepth 태스크 실행 실패
- 파일 경로 불일치

### 문제 4: Python 모듈 에러

```bash
# 필요한 모듈 설치
pip install --user argparse pathlib

# 또는 conda 환경 사용
conda activate hifi-wgs
python3 scripts/generate_qc_report.py
```

---

## 💡 Best Practices

### 1. 정기적인 QC 리포트 생성

배치 작업마다 타임스탬프가 포함된 리포트를 생성하여 히스토리 유지:

```bash
# 타임스탬프 자동 포함
python3 scripts/generate_qc_report.py \
  --output "QC_Report_$(date +%Y%m%d_%H%M%S).html"
```

### 2. 샘플별 QC 체크리스트

리포트 검토 시 확인사항:
- [ ] Mean coverage ≥ 30x
- [ ] Mapping rate ≥ 95%
- [ ] Mean read length ≥ 10kb
- [ ] Mean read quality ≥ Q20
- [ ] SNP Ti/Tv ratio 2.0-2.2 (정상 범위)
- [ ] Het/Hom ratio 1.5-2.0 (정상 범위)

### 3. 버전 관리

리포트와 함께 파이프라인 버전 정보 저장:

```bash
# 버전 정보 파일 생성
cat > /data_4tb/batch_results/pipeline_version.txt << EOF
Pipeline: HiFi-human-WGS-WDL
Version: v3.1.0
Date: $(date)
Samples: $(ls /data_4tb/batch_results | grep -v logs | wc -l)
EOF

# 리포트 생성 시 포함
python3 scripts/generate_qc_report.py --output "QC_Report_v3.1.0_$(date +%Y%m%d).html"
```

### 4. 리포트 백업

```bash
# 주기적 백업
BACKUP_DIR="/data_4tb/qc_reports_archive"
mkdir -p "${BACKUP_DIR}"

cp /data_4tb/batch_results/QC_Report_*.html "${BACKUP_DIR}/"

# 30일 이상 된 리포트 압축
find "${BACKUP_DIR}" -name "QC_Report_*.html" -mtime +30 -exec gzip {} \;
```

---

## 📚 참고 자료

### 관련 문서
- [BATCH_PROCESSING_GUIDE.md](./BATCH_PROCESSING_GUIDE.md) - 배치 처리 전체 가이드
- [MULTI_BAM_GUIDE.md](./MULTI_BAM_GUIDE.md) - 다중 BAM 파일 처리
- [docs/bam_statistics.md](./docs/bam_statistics.md) - BAM 통계 상세

### QC Metrics 기준
- Coverage: GATK Best Practices (30x for WGS)
- Mapping rate: PacBio HiFi 권장 (>95%)
- Read length: PacBio Revio/Sequel II 평균 (15-20kb)
- Ti/Tv ratio: 1000 Genomes Project (2.0-2.2)

### 외부 도구 문서
- [mosdepth](https://github.com/brentp/mosdepth) - Coverage 분석
- [samtools stats](http://www.htslib.org/doc/samtools-stats.html) - BAM 통계
- [DeepVariant](https://github.com/google/deepvariant) - Variant calling
- [pbsv](https://github.com/PacificBiosciences/pbsv) - Structural variant calling

---

## 🔄 업데이트 히스토리

| 날짜 | 버전 | 변경 내용 |
|------|------|-----------|
| 2026-01-30 | 1.0.0 | QC Report 기능 문서화 초안 작성 |

---

## 📞 문의 및 지원

QC 리포트 관련 문제나 개선 제안이 있으시면:
- 이슈 등록: GitHub Issues
- 이메일: (담당자 이메일)
- 문서 개선: Pull Request 환영

---

**다음 단계:**
- [배치 처리 가이드로 돌아가기](./BATCH_PROCESSING_GUIDE.md)
- [파이프라인 실행하기](./README.md#running-the-workflow)
