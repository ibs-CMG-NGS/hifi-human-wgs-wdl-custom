# Batch Processing - 빠른 시작 가이드

> **💡 완전한 가이드**: [BATCH_PROCESSING_GUIDE.md](./BATCH_PROCESSING_GUIDE.md)  
> **📊 QC 리포트**: [QC_REPORT_GUIDE.md](./QC_REPORT_GUIDE.md)

## 📝 사용 순서

### 1단계: 샘플 정보 작성
`samples.csv` 파일을 편집해서 실제 샘플 정보를 추가하세요:

```bash
vim samples.csv
```

형식:
```
sample_id,sex,bam_files
KTY9537,MALE,/data_4tb/pacbio_rawdata/r84285_20260108_080127/1_A01/hifi_reads/m84285_260108_082608_s1.hifi_reads.bc2016.bam
KTY9538,FEMALE,/data_4tb/pacbio_rawdata/sample2/hifi_reads/sample2.bam
```

**여러 BAM 파일을 하나의 샘플로**: 콜론(`:`)으로 구분
```
KTY9539,MALE,/data_4tb/pacbio_rawdata/sample3/cell1.bam:/data_4tb/pacbio_rawdata/sample3/cell2.bam
```

### 2단계: 스크립트에 실행 권한 부여

```bash
chmod +x create_batch_inputs.sh
chmod +x batch_run_optimized.sh
chmod +x monitor_batch.sh
chmod +x collect_results.sh
```

### 3단계: 입력 JSON 파일 생성

```bash
./create_batch_inputs.sh
```

이 명령은 `samples.csv`를 읽어서 각 샘플마다 `batch_inputs/<sample_id>.inputs.json` 파일을 생성합니다.

### 4단계: 배치 실행

**순차 실행 (권장):**
```bash
export CUDA_VISIBLE_DEVICES=1
./batch_run_optimized.sh sequential
```

**특정 샘플만 실행:**
```bash
./batch_run_optimized.sh sequential KTY9537 KTY9538
```

**병렬 실행 (메모리 충분한 경우):**
```bash
./batch_run_optimized.sh parallel
```

### 5단계: 모니터링

**다른 터미널에서 실시간 모니터링:**
```bash
watch -n 30 ./monitor_batch.sh
```

**특정 샘플 로그 확인:**
```bash
tail -f /data_4tb/hifi-human-wgs-wdl-custom/batch_results/logs/KTY9537.log
```

**GPU 모니터링:**
```bash
watch -n 5 nvidia-smi
```

### 6단계: 결과 수집

모든 샘플 완료 후:
```bash
./collect_results.sh
```

결과 확인:
```bash
ls -lh /data_4tb/hifi-human-wgs-wdl-custom/batch_results/summary/
```

---

## 📂 결과 위치

- **각 샘플 결과**: `/data_4tb/hifi-human-wgs-wdl-custom/batch_results/<sample_id>/out/`
- **VCF 파일**: `<sample_id>/out/phased_small_variant_vcf/`
- **BAM 파일**: `<sample_id>/out/merged_haplotagged_bam/`
- **리포트**: `<sample_id>/out/pharmcat_report_html/`
- **통계**: `<sample_id>/out/stats_file/`

---

## 📊 QC 리포트 확인 ⭐

배치 작업이 완료되면 **자동으로 QC 리포트**가 생성됩니다!

```bash
# 리포트 위치 확인
ls -lh /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_*.html

# 브라우저로 열기
firefox /data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report_*.html

# Windows에서 접근
# \\wsl.localhost\Ubuntu\data_4tb\hifi-human-wgs-wdl-custom\batch_results\QC_Report_*.html
```

**리포트 내용:**
- ✅ 전체 샘플 요약 통계
- ✅ Coverage, mapping rate, variant counts
- ✅ QC Pass/Fail 자동 판정

**상세 가이드**: [QC_REPORT_GUIDE.md](./QC_REPORT_GUIDE.md)

---

## ⚠️ 주의사항

1. **디스크 공간**: 샘플당 약 500GB 필요
   ```bash
   df -h /data_4tb
   ```

2. **GPU 설정**: GPU 1번만 사용
   ```bash
   export CUDA_VISIBLE_DEVICES=1
   ```

3. **Conda 환경**: hifi-human-wgs 활성화 확인
   ```bash
   conda activate hifi-human-wgs
   ```

4. **메모리**: 순차 실행 권장 (샘플당 128GB 필요)

---

## 🔧 문제 해결

**중간 파일 정리 (디스크 공간 부족 시):**
```bash
rm -rf /data_4tb/hifi-human-wgs-wdl-custom/batch_results/*/call-*/_miniwdl_*
```

**샘플 하나만 재실행:**
```bash
./batch_run_optimized.sh sequential KTY9538
```

**에러 확인:**
```bash
grep -i "error\|failed" /data_4tb/hifi-human-wgs-wdl-custom/batch_results/logs/*.log
```

---

## 📊 예상 실행 시간

- **단일 샘플 (GPU 모드)**: 7-10시간
- **3개 샘플 (순차)**: 21-30시간
- **재실행 (Call Cache 사용)**: 즉시 완료

---

자세한 내용은 `BATCH_PROCESSING_GUIDE.md`를 참조하세요.
