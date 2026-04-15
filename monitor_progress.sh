#!/bin/bash
# Usage: ./monitor_progress.sh <log_file>
LOG=${1:-"BioSample24.run.log"}

if [ ! -f "$LOG" ]; then
    echo "로그 파일 없음: $LOG"
    exit 1
fi

watch -n 10 bash -c "
echo '================================================'
echo '  miniwdl 진행 상황: $(basename $LOG)'
echo '  $(date)'
echo '================================================'
echo ''

# 완료된 태스크
DONE=\$(grep -c 'task done' '$LOG' 2>/dev/null || echo 0)
# 실행 중인 태스크
RUNNING=\$(grep 'task running' '$LOG' 2>/dev/null | tail -5)
# 에러
ERROR=\$(grep -c '\[error\]\|ERROR\|exception' '$LOG' 2>/dev/null || echo 0)

echo \"  완료된 태스크: \$DONE\"
echo \"  에러 발생: \$ERROR\"
echo ''
echo '--- 현재 실행 중 ---'
grep 'task running' '$LOG' 2>/dev/null | awk '{print \"  >\", \$NF}' | tail -3
echo ''
echo '--- 최근 완료 ---'
grep 'task done' '$LOG' 2>/dev/null | tail -5 | awk '{print \"  v\", \$0}'
echo ''
echo '--- 최근 로그 (5줄) ---'
tail -5 '$LOG'
"
