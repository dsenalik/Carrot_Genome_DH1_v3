#!/bin/bash
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob

# Prepare vcf files for Pixy analysis
# pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data
# https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326
# https://github.com/ksamuk/pixy/blob/master/docs/generating_invar/generating_invar.rst

# configuration
pfx=$(echo $0 | sed -e "s/^.*\///" -e "s/\..*$//")
name="kevinJ"
reference="/vcru_share_s/carrot/dcarv3d/data/DCARv3d.fna"
keybase="/vcru_share_t/box-bioinformatics/keyfiles"
invsegments="06.segments.1000000.txt"   # from 06.makeintervals.sub.sh 1000000
varsegments="06.segments.10000000.txt"  # from 06.makeintervals.sub.sh 10000000
bamdir="/vcru_share_t/bam/0082.04"
tmpdir="/dauc3/tmp/23.kevinJ"
varvcfdir="/vcru_share_t/vcf/0082.06/kevinG"
invvcfdir="/vcru_share_t/vcf/0082.$pfx/$name"
filtvarvcfdir="$tmpdir/variantsfiltered"
mkdir -p "$filtvarvcfdir"
filtinvvcfdir="$tmpdir/invsitesfiltered"
mkdir -p "$filtinvvcfdir"
allvcfdir="$tmpdir/Gcombined"
mkdir -p "$allvcfdir"



############################################################
# Entry Point
############################################################
# combine the two VCFs using bcftools concat
for chr in {1..9} ; do
  mrgvarvcf="/vcru_share_t/vcf/0082.06/kevinG/F2/kevinG.DCARv2_Chr$chr.F2.vcf.gz"
  mrginvvcf="$filtinvvcfdir/DCARv3_Chr$chr.vcf.gz"
  mrgallvcf="$allvcfdir/DCARv3_Chr$chr.vcf.gz"
  log="$allvcfdir/DCARv3_Chr$chr.merge.log"
  if [ ! -s "$mrgallvcf" ] ; then
    echo "+++ Merging variant and invariant Chr$chr  $(date)"
    bcftools concat \
      --allow-overlaps \
      "$mrgvarvcf" "$mrginvvcf" \
      -O z \
      -o "$mrgallvcf" \
      &> "$log" &
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done
wait

for chr in {1..9} ; do
  mrgallvcf="$allvcfdir/DCARv3_Chr$chr.vcf.gz"
  log="$allvcfdir/DCARv3_Chr$chr.merge.log"
  if [ ! -s "$mrgallvcf".tbi ] ; then
    echo "+++ Indexing merged variant and invariant Chr$chr  $(date)"
    tabix -p vcf "$mrgallvcf" &
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done
wait



cp -puv --no-preserve=owner "$0" notebook/
echo "$0 Done  $(date)"
touch "ready"
exit 0



#eof
