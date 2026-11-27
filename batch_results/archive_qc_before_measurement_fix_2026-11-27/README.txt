2026-11-27 측정 버그 수정(a44c334) **이전**에 생성된 raw QC 리포트 원본이다.
파일 날짜(mtime)가 원래 생성 시각이다. 폐기하지 말 것 — LIMS 가 이미 소비한
판정이라 감사 추적에 필요하다.

[진짜 원본 5개] 2026-03-30 ~ 2026-04-13 생성
  QC_Report_rawdata_r84285_20260108_20260113.html
  QC_Report_rawdata_r84285_20260219.html
  QC_Report_rawdata_r84285_20260326.html
  QC_Report_rawdata_r84285_20260402.html
  QC_Report_rawdata_r84285_20260408.html

[원본이 아님 1개]
  NOT-AN-ORIGINAL_r84285_20261124_intermediate_2026-11-27T1041.html
  r84285_20261124 run 의 최초본(2026-11-26 생성)은 같은 날 교체되어 남아
  있지 않다. 이 파일은 수정이 이미 반영된 중간본이라 "수정 전" 기록이
  아니다. 혼동을 막기 위해 이름에 표시했다.
  최초본의 값은 docs/QC_REVISION_2026-11-27.md §3.1 표에 기록되어 있다
  (Control RL 23.0-24.9 kb / Control Reads 3,698-4,162 / P1 54.7-61.7%).

[원본 리포트들이 틀렸던 값]
  - Control Read Length : subread 길이 평균(23-49 kb)을 polymerase read
                          length(45-92 kb) 기준 >=50 kb 에 대조해 전원 WARNING
  - Control Read Count  : subread 개수를 ZMW 기준 수 대신 보고
  - P1                  : ccs_report 기반 근사치 (실측과 최대 +16%p 차이)
  - Control Concordance : 미측정
  - 샘플명              : 20260108/20260113 은 8샘플 중 6샘플이 잘못된 웰에
                          귀속되어 있었다 (20260113 은 5웰 전부 한 칸씩 밀림)

대체본은 같은 이름으로 batch_results/ 에 있다. 변경 내역과 샘플별 판정·귀속
전후 대조는 docs/QC_REVISION_2026-11-27.md 를 보라.
