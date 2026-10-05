#!/usr/bin/env bash

# Demonstration only: update the example paths and resource settings before use.
# Example:
#   bash bash/demo_alignment.sh bash/demo_input /path/to/hisat2/index bash/demo_output
#
# Expected input names:
#   sample_A_R1.fastq.gz
#   sample_A_R2.fastq.gz
#

set -euo pipefail

export input_dir="${1:-bash/demo_input}"
export index_prefix="${2:-/path/to/hisat2/index/prefix}"
export output_dir="${3:-bash/demo_output}"
export threads="${THREADS:-4}"
export jobs="${JOBS:-2}"

if [[ "$index_prefix" == "/path/to/hisat2/index/prefix" ]]; then
    printf 'Set the HISAT2 index prefix as the third argument before running.\n' >&2
    exit 2
fi

if [[ ! -d "$input_dir" ]]; then
    printf 'Input directory not found: %s\n' "$input_dir" >&2
    exit 2
fi

mkdir -p -- "$output_dir"

mapfile -d '' -t r1_files < <(
    find "$input_dir" \
        -maxdepth 1 \
        -type f \
        -name '*_R1.fastq.gz' \
        -print0 |
    sort -z
)

if ((${#r1_files[@]} == 0)); then
    printf 'No *_R1.fastq.gz files found in %s\n' "$input_dir" >&2
    exit 1
fi

# Function to align a sample with HISAT2 and process BAM files
align_sample() (
    local r1="$1"
    local r2="${r1%_R1.fastq.gz}_R2.fastq.gz"
    local index_prefix="$2"
    local output_dir="$3"
    local threads="$4"
    local filename="${r1##*/}"
    local sample="${filename%_R1.fastq.gz}"
    local output_bam="${output_dir}/${sample}.bam"
    local tmpdir

    if [[ ! -f "$r1" ]]; then
        printf 'R1 input not found: %s\n' "$r1" >&2
        return 1
    fi
    if [[ ! -f "$r2" ]]; then
        printf 'Matching R2 input not found: %s\n' "$r2" >&2
        return 1
    fi

    tmpdir=$(mktemp -d) || {
        printf 'Could not create temporary directory\n' >&2
        return 1
    }
    trap 'rm -rf -- "$tmpdir"' EXIT

    printf 'Aligning sample %s\n' "$sample"

    # markdup -r removes reads identified as duplicates from the output BAM.
    hisat2 \
        -x "$index_prefix" \
        -1 "$r1" \
        -2 "$r2" \
        --threads "$threads" |
        samtools sort \
            -n \
            -@ "$threads" \
            -T "$tmpdir/name_sort" \
            - |
        samtools fixmate \
            -m \
            - \
            - |
        samtools sort \
            -@ "$threads" \
            -T "$tmpdir/coordinate_sort" \
            - |
        samtools markdup \
            -r \
            - \
            "$output_bam"

    samtools index "$output_bam"
    printf 'Finished sample %s: %s\n' "$sample" "$output_bam"
)

export -f align_sample

## Run it with a loop
for r1 in "${r1_files[@]}"; do
    align_sample "$r1" "$index_prefix" "$output_dir" "$threads"
done


## Run it with parallel
parallel \
    --null \
    --jobs "$jobs" \
    --halt soon,fail=1 \
    bash -c \
      'align_sample "$@"' \
      _ {} \
      "$index_prefix" \
      "$output_dir" \
      "$threads" \
  ::: "${r1_files[@]}"
