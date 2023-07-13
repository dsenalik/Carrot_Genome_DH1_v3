#!/bin/bash
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob

# pixy analysis

# configuration
pfx=$(echo $0 | sed -e "s/^.*\///" -e "s/\..*$//")

populations="$1"
windowsize="$2"
label="$3"
# populations is either 23.suptbl.33.tsv or 23.suptbl.39.tsv

name="kevinJ"
reference="/vcru_share_s/carrot/dcarv3d/data/DCARv3d.fna"
tmpdir="/dauc3/tmp/$pfx.$name"
#filtvarvcfdir="$tmpdir/variantsfiltered"
#filtinvvcfdir="$tmpdir/invsitesfiltered"
allvcfdir="$tmpdir/Gcombined"
outdir="$tmpdir/Gpixy/$label"
mkdir -p "$outdir"
ncpu=32
chunksize=100000  # default is 100000


until [ -f "ready" ] ; do
  sleep 60
done

for chr in {1..9} ; do
  vcf="$allvcfdir/DCARv3_Chr$chr.vcf.gz"
  outfolder="$outdir/Chr$chr"
  log="$outdir/Chr$chr.log"
  pi="$outfolder/pixy_pi.txt"
  if [ ! -s "$pi" ] ; then
    echo "+++ Pixy chromosome $chr  $(date)"
    pixy \
      --stats pi dxy fst \
      --vcf "$vcf" \
      --populations "$populations" \
      --output_folder "$outfolder" \
      --window_size $windowsize \
      --chunk_size $chunksize \
      --n_cores $ncpu \
      --bypass_invariant_check 'yes' \
      &> "$log"
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    echo "+++ Pixy finished chromosome $chr  $(date)"
  else
    echo "+++ Pixy previously run chromosome $chr  $(date)"
  fi
done



cp -puv --no-preserve=owner "$0" notebook/
echo "$0 Done  $(date)"
exit 0



#eof
