#!/usr/bin/env bash
set -euo pipefail

echo "Starting Railway data sync from URLs..."

if ! command -v railway >/dev/null 2>&1; then
  echo "ERROR: railway CLI not found. Install with: npm install -g @railway/cli"
  exit 1
fi

# Default to the latest validated Zenodo record for HumanProof dataset.
ZENODO_RECORD_ID="${ZENODO_RECORD_ID:-18827087}"
ZENODO_API_BASE="https://zenodo.org/api/records/${ZENODO_RECORD_ID}/files"

LOEUF_URL="${LOEUF_URL:-${ZENODO_API_BASE}/LOEUF_scores.csv.gz/content}"
GENEBASS_PLOF_URL="${GENEBASS_PLOF_URL:-${ZENODO_API_BASE}/genebass_pLoF_filtered.pkl/content}"
MODEL_FINAL_URL="${MODEL_FINAL_URL:-${ZENODO_API_BASE}/model_final.json/content}"
GENE_SHAP_DR_URL="${GENE_SHAP_DR_URL:-${ZENODO_API_BASE}/gene_shap_dr.json/content}"
GENE_METADATA_URL="${GENE_METADATA_URL:-${ZENODO_API_BASE}/gene_metadata.csv/content}"
CELLTYPE_METADATA_URL="${CELLTYPE_METADATA_URL:-${ZENODO_API_BASE}/celltype_metadata.csv/content}"
CELLTYPE_FRACTION_URL="${CELLTYPE_FRACTION_URL:-${ZENODO_API_BASE}/celltype_fraction_expressing.csv/content}"
CELLTYPE_LOG1P_URL="${CELLTYPE_LOG1P_URL:-${ZENODO_API_BASE}/celltype_log1p_mean_expression.csv/content}"

# Expected file sizes in bytes (Zenodo record 18827087). Override if needed.
LOEUF_SIZE="${LOEUF_SIZE:-133005}"
GENEBASS_PLOF_SIZE="${GENEBASS_PLOF_SIZE:-1524965298}"
MODEL_FINAL_SIZE="${MODEL_FINAL_SIZE:-1016758}"
GENE_SHAP_DR_SIZE="${GENE_SHAP_DR_SIZE:-476875641}"
GENE_METADATA_SIZE="${GENE_METADATA_SIZE:-1576435}"
CELLTYPE_METADATA_SIZE="${CELLTYPE_METADATA_SIZE:-48019}"
CELLTYPE_FRACTION_SIZE="${CELLTYPE_FRACTION_SIZE:-90046449}"
CELLTYPE_LOG1P_SIZE="${CELLTYPE_LOG1P_SIZE:-104583376}"

# Optional SHA256 checksums. Strongly recommended for large files.
# If unset, the script will skip checksum validation for that file.
LOEUF_SHA256="${LOEUF_SHA256:-}"
GENEBASS_PLOF_SHA256="${GENEBASS_PLOF_SHA256:-}"
MODEL_FINAL_SHA256="${MODEL_FINAL_SHA256:-}"
GENE_SHAP_DR_SHA256="${GENE_SHAP_DR_SHA256:-}"
GENE_METADATA_SHA256="${GENE_METADATA_SHA256:-}"
CELLTYPE_METADATA_SHA256="${CELLTYPE_METADATA_SHA256:-}"
CELLTYPE_FRACTION_SHA256="${CELLTYPE_FRACTION_SHA256:-}"
CELLTYPE_LOG1P_SHA256="${CELLTYPE_LOG1P_SHA256:-}"

railway ssh "mkdir -p /app/data/safety_model_output/dr /app/data/cellxgene"

fetch_to_remote() {
  local url="$1"
  local dest="$2"
  local expected_size="$3"
  local expected_sha="$4"
  local tmp_dest="${dest}.tmp"

  echo
  echo "Downloading: $url"
  echo "  -> $dest"

  local attempt=1
  local max_attempts=5
  while (( attempt <= max_attempts )); do
    echo "  remote download attempt $attempt/$max_attempts"
    if railway ssh "python3 -c \"import urllib.request as u; req=u.Request('$url', headers={'Accept-Encoding':'identity'}); r=u.urlopen(req, timeout=120); open('$tmp_dest','wb').write(r.read())\""; then
      local remote_size
      remote_size=$(railway ssh "wc -c < '$tmp_dest'" | tr -d '[:space:]')

      if [[ -n "$expected_size" && "$remote_size" != "$expected_size" ]]; then
        echo "  ⚠ size mismatch (got $remote_size, expected $expected_size)"
      else
        break
      fi
    fi

    if (( attempt == max_attempts )); then
      echo "ERROR: failed to download $url after $max_attempts attempts"
      exit 1
    fi

    railway ssh "rm -f '$tmp_dest'" || true
    sleep 3
    ((attempt++))
  done

  if [[ -n "$expected_size" ]]; then
    local remote_size
    remote_size=$(railway ssh "wc -c < '$tmp_dest'" | tr -d '[:space:]')
    if [[ "$remote_size" != "$expected_size" ]]; then
      railway ssh "rm -f '$tmp_dest'"
      echo "ERROR: size mismatch for $dest"
      echo "  expected: $expected_size"
      echo "  actual:   $remote_size"
      exit 1
    fi
    echo "  ✓ size verified"
  fi

  if [[ -n "$expected_sha" ]]; then
    local remote_sha
    remote_sha=$(railway ssh "sha256sum '$tmp_dest' | awk '{print \$1}'" | tr -d '[:space:]')
    if [[ "$remote_sha" != "$expected_sha" ]]; then
      railway ssh "rm -f '$tmp_dest'"
      echo "ERROR: checksum mismatch for $dest"
      echo "  expected: $expected_sha"
      echo "  actual:   $remote_sha"
      exit 1
    fi
    echo "  ✓ checksum verified"
  else
    echo "  ⚠ checksum not provided; integrity verification skipped"
  fi

  railway ssh "mv '$tmp_dest' '$dest'"
  railway ssh "ls -lh '$dest'"
}

fetch_to_remote "$LOEUF_URL" "/app/data/LOEUF_scores.csv.gz" "$LOEUF_SIZE" "$LOEUF_SHA256"
fetch_to_remote "$GENEBASS_PLOF_URL" "/app/data/genebass_pLoF_filtered.pkl" "$GENEBASS_PLOF_SIZE" "$GENEBASS_PLOF_SHA256"

if [[ -n "${MODEL_FINAL_URL:-}" ]]; then
  fetch_to_remote "$MODEL_FINAL_URL" "/app/data/safety_model_output/dr/model_final.json" "$MODEL_FINAL_SIZE" "$MODEL_FINAL_SHA256"
else
  echo
  echo "Skipping model_final.json (MODEL_FINAL_URL not provided)"
fi

fetch_to_remote "$GENE_SHAP_DR_URL" "/app/data/safety_model_output/dr/gene_shap_dr.json" "$GENE_SHAP_DR_SIZE" "$GENE_SHAP_DR_SHA256"
fetch_to_remote "$GENE_METADATA_URL" "/app/data/cellxgene/gene_metadata.csv" "$GENE_METADATA_SIZE" "$GENE_METADATA_SHA256"
fetch_to_remote "$CELLTYPE_METADATA_URL" "/app/data/cellxgene/celltype_metadata.csv" "$CELLTYPE_METADATA_SIZE" "$CELLTYPE_METADATA_SHA256"
fetch_to_remote "$CELLTYPE_FRACTION_URL" "/app/data/cellxgene/celltype_fraction_expressing.csv" "$CELLTYPE_FRACTION_SIZE" "$CELLTYPE_FRACTION_SHA256"
fetch_to_remote "$CELLTYPE_LOG1P_URL" "/app/data/cellxgene/celltype_log1p_mean_expression.csv" "$CELLTYPE_LOG1P_SIZE" "$CELLTYPE_LOG1P_SHA256"

echo
echo "Remote data sync complete."
railway ssh "rm -f /app/data/_stdin_test.txt /app/data/test.txt /app/data/cellxgene/*.tmp /app/data/safety_model_output/dr/*.tmp"
railway ssh "ls -lh /app/data /app/data/cellxgene /app/data/safety_model_output/dr"
railway ssh "wc -c /app/data/LOEUF_scores.csv.gz /app/data/genebass_pLoF_filtered.pkl /app/data/safety_model_output/dr/model_final.json /app/data/safety_model_output/dr/gene_shap_dr.json /app/data/cellxgene/gene_metadata.csv /app/data/cellxgene/celltype_metadata.csv /app/data/cellxgene/celltype_fraction_expressing.csv /app/data/cellxgene/celltype_log1p_mean_expression.csv"
