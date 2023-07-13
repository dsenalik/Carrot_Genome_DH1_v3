#!/bin/bash
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob

# Find the centromere motif in the V3 genome assembly

# configuration
pfx=$(echo $0 | sed -e "s/^.*\///" -e "s/\..*$//")
# page 44 Nov. 26, 2014
marinamotif="CCAACTTCAAACGAGTCAAGAATGAAGCTACAAGTTGTTT"
# page 44 Nov. 28, 2014
consensusmotif="TCAACTTCAAACGAGTCTGGAATGAAGCTACAAGTTGTTT
CCAACTTCAATCTAGYCAATAATGAAGCWACAAGTTGTTT
CCRGGCTCAAACGAGTCAAGAATGAAACTACAAGTTGATT
CYAATTTYAAAC AACCGAAAATGAAGCTACAAGTTGTTT"
# master assembly
mastergff='/vcru_share_s/carrot/dcarv3d/data/DCARv3d.master.gff3'
# chromosome lengths (from $mastergff)
chrlen=(0 56320272 59175025 64558952 51152908 47556597 40885772 41844801 33368437 45850561)
# configuration of bin size for this analysis
binsize=250000

# input files
v3fasta="/vcru_share_s/carrot/dcarv3d/data/DCARv3d.fna"
v3db="/vcru_share_s/blastdbs/DCARv3d"

# output files
querymarina="$pfx.marinamotif.fna"
outmarina="$pfx.marinamotif.tsv"
queryconsensus="$pfx.consensusmotif.fna"
outconsensus="$pfx.consensusmotif.tsv"


# store two copies of the motif sequences so that there is no edge effect

if [ ! -s "$querymarina" ] ; then
  echo "Creating motif fasta \"$querymarina\"  $(date)"
  echo ">cent1
$marinamotif
$marinamotif" > "$querymarina"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
fi

if [ ! -s "$queryconsensus" ] ; then
  echo "Creating motif fasta \"$queryconsensus\"  $(date)"
  echo ">cent1
$consensusmotif
$consensusmotif" > "$queryconsensus"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
fi



if [ ! -s "$outmarina" ] ; then
  echo "Blast Marina  $(date)"
  blastn \
    -db "$v3db" \
    -query "$querymarina" \
    -out "$outmarina" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp" \
    -dust no \
    -max_target_seqs 50000
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
fi
  


if [ ! -s "$outconsensus" ] ; then
  echo "Blast consensus  $(date)"
  blastn \
    -db "$v3db" \
    -query "$queryconsensus" \
    -out "$outconsensus" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp" \
    -dust no \
    -max_target_seqs 50000
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
fi



for chr in {1..9} ; do
  echo "Summarizing chromosome $chr  $(date)"
  # Length of current chromosome, defined in config section
  len=${chrlen[$chr]}

  # Create file of bins based on chromosome length
  perl - $chr $len $binsize <<'  __EOP__'
  use strict;
  use warnings;
  my $chr = shift;
  my $len = shift;
  my $binsize = shift; # defined in parent bash script config section
  my $i = 0;
  open (my $OUTF, '>', 'bins.bed');
  while ($i < $len)
    {
      my $e = $i + $binsize;
      if ($e > $len) { $e = $len; }
      print $OUTF join("\t", "DCARv3_Chr$chr", $i, $e), "\n";
      $i += $binsize;
    }
  __EOP__
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi

  # Extract blast hit data just for this chromosome, convert to bed coordinate system
  cat "$outconsensus" \
    | perl -lanF/"\t"/ -e 'if ($F[8]>$F[9]){($F[8],$F[9])=($F[9],$F[8])} print join ("\t", $F[1], $F[8]-1, $F[9])' \
    | grep "^DCARv3_Chr$chr" \
    | sort-bed - \
    > "$pfx.tmp.1"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi

  bedtools coverage -a bins.bed -b "$pfx.tmp.1" > "$pfx.DCARv3_Chr$chr.tsv"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi

  # preprocess for genome plotting. Convert fraction coverage (0-1) to percent (0-100)
  echo "Chr	Pos	Value" > "$pfx.tmp.2"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  cat "$pfx.DCARv3_Chr$chr.tsv" \
    | perl -lanF/"\t"/ -e 'print join ("\t", $F[0], $F[1], $F[6]*100)' \
    >> "$pfx.tmp.2"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  bb.genomeplotter \
    --infile="$pfx.tmp.2" \
    --outfile="notebook/$pfx.DCARv3_Chr$chr.png" \
    --mastergff="$mastergff" \
    --chr="DCARv3_Chr$chr" \
    --bars --filled \
    --legend="DCARv3 Chr$chr" \
    --nolabels \
    --binsize=$binsize \
    --ymax=100
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  rm -f bins.bed $pfx.tmp.1 $pfx.tmp.2
done




cp -puv --no-preserve=owner "$0" notebook/
echo "$0 Done  $(date)"
exit 0



#eof
