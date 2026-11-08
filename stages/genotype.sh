#!/usr/bin/env bash
# Один GenotypeGVCFs на вид, все выбранные аутосомы в одном вызове.
# Usage: genotype.sh REF.fa GENOMICS_DB_DIR OUT.bcf AUTOSOME [AUTOSOME ...]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/resource_lib.bash"
if [[ "${1:-}" == '-h' || "${1:-}" == '--help' ]]; then
    echo 'Usage: genotype.sh REF.fa GENOMICS_DB_DIR OUT.bcf AUTOSOME [AUTOSOME ...]'; exit 0
fi
if (( $# < 4 )); then echo 'Usage: genotype.sh REF.fa GENOMICS_DB_DIR OUT.bcf AUTOSOME [AUTOSOME ...]' >&2; exit 2; fi
reference_fasta=$1
db_dir=$2
output_bcf=$3
shift 3
[[ -f "$db_dir/.completed" ]] || { echo "Не завершён GenomicsDB: $db_dir" >&2; exit 2; }
[[ ! -e "$output_bcf" && ! -e "${output_bcf}.csi" ]] || {
    echo "Выходной BCF уже существует: $output_bcf" >&2; exit 2;
}
threads=${MIN_SNP_THREADS:-1}
minsnp_check_threads "$threads"
heap_mb=$(minsnp_heap_mb 75 16384)
gc_threads=$((threads < 8 ? threads : 8))
mkdir -p "$(dirname "$output_bcf")"
intervals=()
for contig in "$@"; do intervals+=(-L "$contig"); done
# GATK пишет VCF.gz, затем преобразуем его в BCF. Не оставляем временный VCF.
tmp_vcf="${output_bcf}.incomplete.vcf.gz"
tmp_bcf="${output_bcf}.incomplete.bcf"
trap 'rm -f -- "$tmp_vcf" "${tmp_vcf}.tbi" "${tmp_vcf}.csi" "$tmp_bcf" "${tmp_bcf}.csi"' EXIT
[[ ! -e "$tmp_vcf" && ! -e "$tmp_bcf" ]] || {
    echo 'Остался неполный файл от старого GenotypeGVCFs; проверьте его перед повторным запуском' >&2; exit 2;
}
gatk --java-options "-Xmx${heap_mb}m -XX:ActiveProcessorCount=${threads} -XX:ParallelGCThreads=${gc_threads}" \
    GenotypeGVCFs -R "$reference_fasta" -V "gendb://${db_dir}" \
    "${intervals[@]}" -O "$tmp_vcf"
test -s "$tmp_vcf"
bcftools view -Ob -o "$tmp_bcf" "$tmp_vcf"
bcftools index -c "$tmp_bcf"
bcftools index -n "$tmp_bcf" >/dev/null
mv -- "$tmp_bcf" "$output_bcf"
mv -- "${tmp_bcf}.csi" "${output_bcf}.csi"
trap - EXIT
rm -f -- "$tmp_vcf" "${tmp_vcf}.tbi" "${tmp_vcf}.csi"
