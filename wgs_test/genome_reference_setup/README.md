# Reference preparation — historical helper

This directory retains `setup_genomes.sh` and `reference_summary.py`, originally written for the earlier BWA-MEM workflow. **Do not re-run reference downloads or indexing for `pilot_20261030`.** The current references and BWA-MEM3 indices were prepared and successfully used on `sphinx` before the 02.11.2026 snapshot.

| Species | Assembly | Autosomes |
|---|---|---:|
| *Delphinapterus leucas* | `GCA_059928025.1` (BTC_Dleu_Yukka2) | 21 |
| *Ursus maritimus* | `GCA_062004589.1` (BTC_Umar_Aika3) | 36 |

**Current source of truth:** `config/references.yaml` (absolute FASTA paths and exact contig list). Required BWA-MEM3 index components include `<FASTA>.amb`, `.ann`, `.pac` and **`.bwt.2bit.64`**. The historic helper's `bwa index` outputs (`.bwt`, `.sa`) are **not** substitutes for `bwa-mem3 index`.

The cluster contains nodes lacking AVX2, required by installed BWA-MEM3. The pilot's reference indexing was moved to the AVX2-capable `meduza-1`. Any future re-indexing requires an explicit independent task with sufficient RAM, a verified FASTA and an appropriate compute node, not a command run automatically by a setup guide.

See [technical guide](../../docs/PROJECT_GUIDE.md) for current environment, reference paths and non-destructive checks. This helper is retained for provenance; it should be audited separately before use against a new assembly.
