#!/usr/bin/env bash

# Java heap: доля выделенной памяти с верхним ограничением.
minsnp_heap_mb() {
    local heap=$(( ${MIN_SNP_MEMORY_MB:?} * $1 / 100 ))
    printf '%s\n' "$(( heap < $2 ? heap : $2 ))"
}
