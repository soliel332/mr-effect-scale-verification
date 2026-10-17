#!/usr/bin/env bash
set -euo pipefail

analysis_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest_file="$analysis_root/ldsc_input_manifest.tsv"
download_dir="$analysis_root/private_raw_inputs"

mkdir -p "$download_dir"

download_one() {
  analysis_label="$1"
  accession="$2"
  source_url="$3"
  local_filename="$4"
  destination_file="$download_dir/$local_filename"
  echo "DOWNLOAD_START $analysis_label $accession $local_filename"
  /usr/bin/curl \
    --fail \
    --location \
    --retry 5 \
    --retry-delay 5 \
    --retry-all-errors \
    --continue-at - \
    --output "$destination_file" \
    "$source_url"
  echo "DOWNLOAD_DONE $analysis_label $accession $local_filename"
}

# Keep downloads sequential. This makes a failed transfer stop the script
# immediately and prevents the parent shell from reporting completion before a
# background job launched inside a pipeline subshell has finished.
while IFS=$'\t' read -r analysis_label manuscript_role accession trait source_url local_filename reported_n ancestry build distribution_note; do
  download_one "$analysis_label" "$accession" "$source_url" "$local_filename"
done < <(tail -n +2 "$manifest_file")

shasum -a 256 "$download_dir"/* > "$analysis_root/private_raw_inputs_SHA256SUMS.txt"
echo "ALL_DOWNLOADS_DONE"
