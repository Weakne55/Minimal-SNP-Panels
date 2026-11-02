#!/usr/bin/env bash
# One-time reference installation for beluga and polar bear WGS.
# No reference generation/indexing is performed by the WGS pipeline itself.
set -euo pipefail

BELUGA_ACC='GCA_059928025.1'
POLAR_ACC='GCA_062004589.1'
BELUGA_DIR='Delphinapterus_leucas_GCA_059928025.1'
POLAR_DIR='Ursus_maritimus_GCA_062004589.1'

usage() {
  cat <<'HELP'
Usage:
  bash setup_genomes.sh download --root /mnt/tank/scratch/plem/data/genomes
  bash setup_genomes.sh index    --root /mnt/tank/scratch/plem/data/genomes
  bash setup_genomes.sh check    --root /mnt/tank/scratch/plem/data/genomes
  bash setup_genomes.sh inspect  --root /mnt/tank/scratch/plem/data/genomes
  bash setup_genomes.sh all      --root /mnt/tank/scratch/plem/data/genomes
  bash setup_genomes.sh index-one /absolute/path/to/genome.fa

Download uses NCBI Datasets for BOTH assemblies, as requested by the team.
No independent Dryad reference, no extra haplotypes and no annotation files.
Index creates only missing .fai, .dict and BWA files, once, on the server.
The WGS pipeline itself is never modified or run by this script.
HELP
}

err() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || err "Missing $1 in PATH"; }
msg() { printf '[reference-setup] %s\n' "$*" >&2; }

[[ $# -ge 1 ]] || { usage; exit 2; }
action=$1; shift
if [[ $action == index-one ]]; then
  [[ $# -eq 1 ]] || { usage; exit 2; }
  single_fasta=$1
else
  root=''
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --root) [[ $# -ge 2 ]] || err '--root requires a directory'; root=$2; shift 2;;
      -h|--help) usage; exit 0;;
      *) err "Unknown option: $1";;
    esac
  done
  [[ -n $root ]] || err 'Specify --root for shared genomes directory'
  mkdir -p -- "$root"
  root=$(cd "$root" && pwd -P)
fi

check_fasta() {
  [[ -s $1 ]] || err "FASTA does not exist or is empty: $1"
  grep -qm 1 '^>' "$1" || err "FASTA header not found: $1"
}

beluga_fasta() {
  local folder="$root/$BELUGA_DIR/ncbi_dataset/data/$BELUGA_ACC"
  [[ -d $folder ]] || return 0
  find "$folder" -maxdepth 1 -type f -name '*_genomic.fna' -print -quit
}
polar_fasta() {
  local folder="$root/$POLAR_DIR/ncbi_dataset/data/$POLAR_ACC"
  [[ -d $folder ]] || return 0
  find "$folder" -maxdepth 1 -type f -name '*_genomic.fna' -print -quit
}

record() {
  local fasta=$1 source=$2 folder
  folder=$(dirname "$fasta")
  printf '%s\n' "$source" > "$folder/reference_source.url.txt"
  (cd "$folder" && sha256sum "$(basename "$fasta")" > "$(basename "$fasta").sha256")
  msg "SHA256 saved: $fasta.sha256"
}

download_one() {
  local species=$1 accession=$2 directory=$3 fasta_fn=$4
  local base="$root/$directory" zip fasta
  mkdir -p "$base"
  fasta=$($fasta_fn)
  if [[ -n $fasta && -s $fasta ]]; then
    msg "$species FASTA already present: $fasta"
    return
  fi
  need datasets; need unzip; need sha256sum
  zip="$base/ncbi_dataset.zip"
  if [[ ! -s $zip ]]; then
    msg "Downloading $species reference $accession from NCBI..."
    datasets download genome accession "$accession" --include genome,seq-report \
      --filename "$zip"
  else
    msg "Reusing existing NCBI dataset ZIP: $zip"
  fi
  unzip -tq "$zip" >/dev/null || err "Corrupt NCBI ZIP: $zip (remove it and retry)"
  unzip -n -q "$zip" -d "$base"
  fasta=$($fasta_fn)
  [[ -n $fasta ]] || err "NCBI package does not contain $accession/*_genomic.fna"
  check_fasta "$fasta"
  record "$fasta" "NCBI Assembly https://www.ncbi.nlm.nih.gov/datasets/genome/$accession/"
  msg "$species FASTA: $fasta"
  msg "Original archive kept at $zip (can remove after verification to save space)"
}

download_beluga() { download_one Beluga "$BELUGA_ACC" "$BELUGA_DIR" beluga_fasta; }
download_polar() { download_one 'Polar bear' "$POLAR_ACC" "$POLAR_DIR" polar_fasta; }

reference_dictionary() {
  local fa=$1 stem
  stem=${fa%.*}
  printf '%s.dict\n' "$stem"
}

index_one() {
  local fa=$1 dict x need_bwa=0
  check_fasta "$fa"
  need samtools; need bwa; need gatk
  dict=$(reference_dictionary "$fa")
  if [[ ! -s "$fa.fai" ]]; then
    msg "Creating FASTA index: $fa.fai"
    samtools faidx "$fa"
  else msg "Existing FASTA index: $fa.fai"; fi
  if [[ ! -s $dict ]]; then
    msg "Creating GATK dictionary: $dict"
    gatk CreateSequenceDictionary -R "$fa" -O "$dict"
  else msg "Existing GATK dictionary: $dict"; fi
  for x in amb ann bwt pac sa; do
    [[ -s "$fa.$x" ]] || need_bwa=1
  done
  if (( need_bwa )); then
    msg "Creating BWA-MEM index for $fa (one-time operation)..."
    bwa index -a bwtsw "$fa"
  else msg "Existing BWA-MEM index: $fa.{amb,ann,bwt,pac,sa}"; fi
  check_one "$fa"
}

check_one() {
  local fa=$1 dict x
  check_fasta "$fa"
  dict=$(reference_dictionary "$fa")
  for x in "$fa.fai" "$dict" "$fa.amb" "$fa.ann" "$fa.bwt" "$fa.pac" "$fa.sa"; do
    [[ -s $x ]] || err "Missing or empty required index: $x"
  done
  python3 "$(dirname "${BASH_SOURCE[0]}")/reference_summary.py" check "$fa"
}

get_pair() {
  local b p
  b=$(beluga_fasta); p=$(polar_fasta)
  [[ -n $b && -s $b ]] || err "Beluga FASTA missing; run download first ($BELUGA_ACC)"
  [[ -n $p && -s $p ]] || err "Polar bear FASTA missing ($POLAR_ACC); run download first"
  printf '%s\n%s\n' "$b" "$p"
}

case "$action" in
  download) download_beluga; download_polar;;
  index)
    mapfile -t fasta_pair < <(get_pair)
    [[ ${#fasta_pair[@]} -eq 2 ]] || err 'Both assemblies must be downloaded before indexing/checking'
    for fa in "${fasta_pair[@]}"; do index_one "$fa"; done;;
  check)
    mapfile -t fasta_pair < <(get_pair)
    [[ ${#fasta_pair[@]} -eq 2 ]] || err 'Both assemblies must be downloaded before indexing/checking'
    for fa in "${fasta_pair[@]}"; do check_one "$fa"; done;;
  inspect)
    mapfile -t fasta_pair < <(get_pair)
    [[ ${#fasta_pair[@]} -eq 2 ]] || err 'Both assemblies must be downloaded before indexing/checking'
    for fa in "${fasta_pair[@]}"; do
      python3 "$(dirname "${BASH_SOURCE[0]}")/reference_summary.py" inspect "$fa" --limit 45
    done;;
  all)
    download_beluga; download_polar
    mapfile -t fasta_pair < <(get_pair)
    [[ ${#fasta_pair[@]} -eq 2 ]] || err 'Both assemblies must be downloaded before indexing/checking'
    for fa in "${fasta_pair[@]}"; do index_one "$fa"; done
    for fa in "${fasta_pair[@]}"; do
      python3 "$(dirname "${BASH_SOURCE[0]}")/reference_summary.py" inspect "$fa" --limit 45
    done;;
  index-one) index_one "$single_fasta";;
  -h|--help|help) usage;;
  *) usage; exit 2;;
esac
