"""Единый бюджет CPU и оперативной памяти для локального запуска/SLURM."""
from __future__ import annotations

import os
import re


def memory_mib(value: str) -> int:
    """Преобразовать SLURM-память `64G`, `32768M` или `65536` в MiB."""
    match = re.fullmatch(r"([1-9][0-9]*)([MGT]?)", str(value).strip().upper())
    if not match:
        raise ValueError(f"Неверное значение памяти: {value!r} (пример: 64G)")
    number = int(match.group(1))
    factor = {"": 1, "M": 1, "G": 1024, "T": 1024 * 1024}[match.group(2)]
    return number * factor


def configure_resources(args, cfg: dict) -> dict:
    """Проверить ресурсы и учесть ограничения текущего SLURM allocation.

    Используем реально выделенные CPU, даже если config.yaml требует больше.
    Не создаём пул процессов и не предполагаем, что RAM линейна числу CPU.
    """
    allocation = int(os.environ.get("SLURM_CPUS_PER_TASK", "0") or 0)
    requested_threads = (args.threads if args.threads is not None
                         else allocation if allocation else cfg.get("threads", 1))
    if type(requested_threads) is not int or requested_threads < 1:
        raise ValueError("--threads / threads в YAML должен быть целым числом >= 1")
    requested_memory = str(args.memory or cfg["slurm"]["memory"]).upper()
    available_memory = memory_mib(requested_memory)
    if available_memory < 1024:
        raise ValueError("Для работы требуется не менее 1024 MiB памяти")

    allocated_threads = allocation
    threads = min(requested_threads, allocated_threads) if allocated_threads else requested_threads
    allocated_memory = int(os.environ.get("SLURM_MEM_PER_NODE", "0") or 0)
    if not allocated_memory and os.environ.get("SLURM_MEM_PER_CPU"):
        allocated_memory = int(os.environ["SLURM_MEM_PER_CPU"]) * threads
    if allocated_memory:
        available_memory = min(available_memory, allocated_memory)
    if available_memory < 1024:
        raise ValueError("В выделенном SLURM job менее 1024 MiB RAM")

    cfg["threads"] = threads
    cfg["_memory_mb"] = available_memory
    cfg["_requested_threads"] = requested_threads
    cfg["_requested_memory"] = requested_memory
    os.environ["MIN_SNP_THREADS"] = str(threads)
    os.environ["MIN_SNP_MEMORY_MB"] = str(available_memory)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    print(f"Ресурсы: CPU={threads} RAM={available_memory} MiB (запрос: {requested_threads}, {requested_memory})")
    return cfg
