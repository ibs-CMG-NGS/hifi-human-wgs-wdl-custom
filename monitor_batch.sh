#!/bin/bash
# monitor_batch.sh
# Batch 처리 중 리소스 모니터링

BATCH_RESULTS="/data_4tb/hifi-human-wgs-wdl-custom/batch_results"
LOG_DIR="${BATCH_RESULTS}/logs"

echo "========================================"
echo "Batch Processing Monitor"
echo "========================================"
echo ""

# 실행 중인 miniwdl 프로세스
echo "Running processes:"
ps aux | grep -E "(miniwdl|pbmm2|deepvariant)" | grep -v grep || echo "  None"
echo ""

# CPU 및 메모리 사용량
echo "System resources:"
echo "  CPU: $(top -bn1 | grep "Cpu(s)" | awk '{print $2}')% used"
echo "  Memory:"
free -h | grep -E "Mem|Swap"
echo ""

# GPU 상태
echo "GPU status:"
nvidia-smi --query-gpu=index,name,temperature.gpu,utilization.gpu,memory.used,memory.total --format=csv,noheader,nounits 2>/dev/null || echo "  Not available"
echo ""

# 디스크 사용량
echo "Disk usage:"
df -h / /data_4tb
echo ""

# 샘플별 상태
echo "Sample status:"
for log in ${LOG_DIR}/*.log; do
    if [[ -f "$log" ]]; then
        sample=$(basename "$log" .log)
        
        # 마지막 로그 라인
        last_line=$(tail -1 "$log" 2>/dev/null)
        
        # 완료 여부 확인
        if grep -q "done" "$log" 2>/dev/null; then
            status="✓ COMPLETED"
        elif grep -q "error\|failed" "$log" 2>/dev/null; then
            status="✗ FAILED"
        elif [[ -f "${LOG_DIR}/${sample}.pid" ]]; then
            pid=$(cat "${LOG_DIR}/${sample}.pid")
            if ps -p $pid > /dev/null 2>&1; then
                status="⏳ RUNNING (PID: $pid)"
            else
                status="⏹ STOPPED"
            fi
        else
            status="⏸ WAITING"
        fi
        
        echo "  ${sample}: ${status}"
    fi
done
echo ""
echo "========================================"
