version 1.0

import "wdl-common/wdl/structs.wdl"
import "wdl-common/wdl/tasks/base_modification_per_sample.wdl" as BMS

# PacBio HiFi per-sample 5hmC/5mC 경로 (replicate 단위 통계 분석용)
#
# base_modification_known.wdl은 그룹 내 BAM을 병합(pooling)해 replicate 정보가
# 사라지므로, 샘플 하나씩 처리해 DSS 등 replicate 기반 통계의 입력을 만든다.
#
# 흐름: raw HiFi BAM → jasmine(미정렬, kinetics 제거) → MM 필터 + pbmm2 정렬
#       → modkit CpG pileup(공통 고정 threshold) → 5mC/5hmC per-sample TSV
#
# kinetics 보존 정렬(원본×1.28)을 거치지 않아 디스크 사용량이 크게 줄어든다.
# 이미 jasmine 호출·정렬된 BAM이 있으면(예: 파일럿 1차 run) existing_modbam으로
# 넣어 pileup만 수행한다.
#
# 6mA는 WGS BAM의 장비 호출값을 쓰므로 이 워크플로가 아니라
# base_mod_per_sample/run_6mA_pileup.sh에서 처리한다.

workflow base_modification_per_sample {
  meta {
    description: "HiFi per-sample 5mC/5hmC: jasmine(미정렬 raw) → pbmm2 → modkit CpG pileup"
  }

  input {
    String sample_id

    File?  raw_hifi_bam
    File?  existing_modbam
    File?  existing_modbai

    File   ref_fasta
    File   ref_fai

    Int    min_coverage     = 5
    Float  filter_threshold = 0.7

    RuntimeAttributes default_runtime_attributes
  }

  if (!defined(existing_modbam)) {
    call BMS.jasmine_unaligned {
      input:
        raw_hifi_bam       = select_first([raw_hifi_bam]),
        sample_id          = sample_id,
        runtime_attributes = default_runtime_attributes
    }

    call BMS.pbmm2_align_modbam {
      input:
        jasmine_bam        = jasmine_unaligned.jasmine_bam,
        ref_fasta          = ref_fasta,
        ref_fai            = ref_fai,
        sample_id          = sample_id,
        runtime_attributes = default_runtime_attributes
    }
  }

  File modbam = select_first([pbmm2_align_modbam.modbam, existing_modbam])
  File modbai = select_first([pbmm2_align_modbam.modbai, existing_modbai])

  call BMS.modbam_cpg_pileup {
    input:
      modbam             = modbam,
      modbai             = modbai,
      ref_fasta          = ref_fasta,
      ref_fai            = ref_fai,
      sample_id          = sample_id,
      min_coverage       = min_coverage,
      filter_threshold   = filter_threshold,
      runtime_attributes = default_runtime_attributes
  }

  output {
    File? jasmine_bam    = jasmine_unaligned.jasmine_bam
    File? aligned_modbam = pbmm2_align_modbam.modbam
    File? aligned_modbai = pbmm2_align_modbam.modbai
    File  fivemC_tsv     = modbam_cpg_pileup.fivemC_tsv
    File  fivehmC_tsv    = modbam_cpg_pileup.fivehmC_tsv
    File  modkit_log     = modbam_cpg_pileup.modkit_log
  }
}
