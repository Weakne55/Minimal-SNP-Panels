#!/usr/bin/env bash
# Общие функции CPU/RAM для стадий FASTQ -> VCF.
# Источник памяти: MIN_SNP_MEMORY_MB (MiB), установленный Python CLI.

minsnp_check_threads() {
    local threads=$1
    if [[ ! "$threads" =~ ^[0-9]+$ ]] || (( threads < 1 )); then
        echo "ОШИБКА: THREADS должен быть целым числом >= 1" >&2
        exit 2
    fi
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]] && (( threads > SLURM_CPUS_PER_TASK )); then
        echo "ОШИБКА: запрошено $threads CPU, SLURM выделил ${SLURM_CPUS_PER_TASK}" >&2
        exit 2
    fi
}

minsnp_memory_mb() {
    local memory=${MIN_SNP_MEMORY_MB:-65536}
    if [[ ! "$memory" =~ ^[0-9]+$ ]] || (( memory < 1024 )); then
        echo "ОШИБКА: MIN_SNP_MEMORY_MB должен быть целым числом >= 1024" >&2
        exit 2
    fi
    printf '%s\n' "$memory"
}

# Один основной поток в samtools плюс указанное число дополнительных (-@).
minsnp_samtools_workers() {
    local threads=$1
    local workers=$((threads > 1 ? threads - 1 : 0))
    (( workers > 32 )) && workers=32
    printf '%s\n' "$workers"
}

# Java heap с запасом под GC, native memory, файловые буферы.
# arg1: доля от общего лимита (в процентах), arg2: верхняя граница MiB.
minsnp_heap_mb() {
    local percent=$1 maximum=$2
    local memory heap
    memory=$(minsnp_memory_mb)
    heap=$((memory * percent / 100))
    (( heap > maximum )) && heap=$maximum
    (( heap < 512 )) && heap=512
    printf '%s\n' "$heap"
}
