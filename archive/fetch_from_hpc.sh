#!/bin/bash
# fetch_from_hpc.sh  --  RUN THIS ON THE HPC (ip-172-25-88-203)
# Step 1 measures, Step 2 transfers. Nothing moves until you run step 2.

ROOT=/mnt/lustre/lab/gwang/mmcloud_2026/ftp_fgc_xqtl
DEST=<user>@scc1.bu.edu:${AD_LOCI_ROOT:?set AD_LOCI_ROOT}/
LIST=${1:-/tmp/transfer_list_minimal.txt}

cd "$ROOT" || { echo "cannot cd $ROOT"; exit 1; }
: > /tmp/xfer_ok.txt; : > /tmp/xfer_missing.txt
total=0
while IFS= read -r p; do
  [ -z "$p" ] && continue
  if [ -e "$p" ]; then
    echo "$p" >> /tmp/xfer_ok.txt
    sz=$(du -sb "$p" 2>/dev/null | cut -f1)
    total=$(( total + ${sz:-0} ))
  else
    echo "$p" >> /tmp/xfer_missing.txt
  fi
done < "$LIST"

echo "=============================================="
echo "list:    $LIST  ($(wc -l < "$LIST") entries)"
echo "found:   $(wc -l < /tmp/xfer_ok.txt)"
echo "missing: $(wc -l < /tmp/xfer_missing.txt)   -> /tmp/xfer_missing.txt"
awk -v b="$total" 'BEGIN{ printf "SIZE:    %.2f GB\n", b/1073741824 }'
echo "=============================================="
echo
echo "10 largest entries:"
while IFS= read -r p; do du -sh "$p" 2>/dev/null; done < /tmp/xfer_ok.txt | sort -rh | head -10
echo
echo "If the size looks acceptable, transfer with:"
echo "  rsync -avh --progress --files-from=/tmp/xfer_ok.txt $ROOT/ $DEST"