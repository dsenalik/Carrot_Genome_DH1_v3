#!/bin/bash
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob

# comparative analysis of DCARv3d with v2 and v3c assemblies

# configuration
pfx=$(echo $0 | sed -e "s/^.*\///" -e "s/\..*$//")
chrs="Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 MT PT"

# input files
v2seq="/vcru_share_s/carrot/LNRQ01/data/130.DHv2.fna"
v3cseq="/vcru_share_s/carrot/dcarv3c/data/DCARv3c.fna.gz"
v3dseq="/vcru_share_s/carrot/dcarv3d/data/DCARv3d.fna.gz"

# output files
outdir="$pfx.outputdir"
nboutdir="notebook/$pfx.outputdir"



# 1. Initialization
for d in $outdir $nboutdir ; do
  if [ ! -d $d ]; then
    echo "1. Creating directory \"$d\"  $(date)"
    mkdir "$d"
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done



# 2. Extract subsequences
for chr in $chrs ; do
  oldchr="/vcru_share_s/carrot/LNRQ01/data/130.DHv2_$chr.fna.gz"
  v3cchr="/vcru_share_s/carrot/dcarv3c/data/DCARv3c_${chr}.fna.gz"
  v3dchr="/vcru_share_s/carrot/dcarv3d/data/DCARv3d_${chr}.fna.gz"

  if [ ! -s "$oldchr" ]; then
    echo "Extracting v2 $chr  $(date)"
    bb.fastagrep --grep="DCARv2_$chr" --infile="$v2seq" --outfile="$oldchr"
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
  if [ ! -s "$v3cchr" ]; then
    echo "Extracting v3c $chr  $(date)"
    bb.fastagrep --grep="DCARv3_$chr" --infile="$v3cseq" --outfile="$v3cchr"
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
  if [ ! -s "$v3dchr" ]; then
    echo "Extracting v3d $chr  $(date)"
    bb.fastagrep --grep="DCARv3_$chr" --infile="$v3dseq" --outfile="$v3dchr"
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done  # for chr



# 3. Comparison with MUMmer v2 chromosome vs. v3d chromosome
for chr in $chrs ; do
  oldchr="/vcru_share_s/carrot/LNRQ01/data/130.DHv2_$chr.fna.gz"
  v3cchr="/vcru_share_s/carrot/dcarv3c/data/DCARv3c_${chr}.fna.gz"
  v3dchr="/vcru_share_s/carrot/dcarv3d/data/DCARv3d_${chr}.fna.gz"
  mummerout3a="$outdir/$pfx.mummer-v2-v3d.$chr"
  mummerout3b="$outdir/$pfx.mummer-v3c-v3d.$chr"

  if [ ! -s "${mummerout3a}1-0.tn.jpg" ]; then
    echo "3a. Comparing v2-v3d $chr with MUMmer  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$oldchr" \
      --outfile="$mummerout3a" \
      --size=large \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --savegnuplot \
      --nopoints \
      --minlength=5000 \
      --minpercent=96 \
      --showgaps=1000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    tn "${mummerout3a}1-0.png" 6
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    cp -puv --no-preserve=owner "${mummerout3a}1-0.png" "${mummerout3a}1-0.tn.jpg" "$nboutdir"/
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi

  if [ ! -s "${mummerout3b}1-0.tn.jpg" ]; then
    echo "3b. Comparing v3c-v3d $chr with MUMmer  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$v3cchr" \
      --outfile="$mummerout3b" \
      --size=large \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --savegnuplot \
      --nopoints \
      --minlength=5000 \
      --minpercent=96 \
      --showgaps=1000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    tn "${mummerout3b}1-0.png" 6
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    cp -puv --no-preserve=owner "${mummerout3b}1-0.png" "${mummerout3b}1-0.tn.jpg" "$nboutdir"/
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done  # for chr
exit



# 4. Comparison with MUMmer organellar with alignment output
# ( organellar are the same for v3c and v3d )
for chr in MT PT ; do
  oldchr="/vcru_share_s/carrot/LNRQ01/data/130.DHv2_$chr.fna.gz"
  v3dchr="/vcru_share_s/carrot/dcarv3d/data/DCARv3d_${chr}.fna.gz"
  mummerout4a="$outdir/$pfx.mummer4a.$chr"

  if [ ! -s "${mummerout4a}1-0.alignment.txt" ]; then
    echo "4a. Comparing $chr with MUMmer  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$oldchr" \
      --outfile="$mummerout4a" \
      --size=none \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --showaligns \
      --minlength=5000 \
      --minpercent=96 \
      --showgaps=1000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done  # for chr



# 5. New region analysis - compare a chromosome to entire old genome = SLOW!
for chr in $chrs ; do
  v3dchr="/vcru_share_s/carrot/dcarv3d/data/DCARv3d_${chr}.fna.gz"
  mummerout5a="$outdir/$pfx.mummer5a-v3d-to-v2.$chr"
  mummerout5b="$outdir/$pfx.mummer5b-v3d-to-v3c.$chr"

  if [ ! -s "${mummerout5a}1-0.tn.jpg" ]; then
    echo "5a. Comparing v3d $chr to v2 genome with MUMmer  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$v2seq" \
      --outfile="$mummerout5a" \
      --size=large \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --savegnuplot \
      --nopoints \
      --minlength=5000 \
      --minpercent=96 \
      --showgaps=100000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    tn "${mummerout5a}1-0.png" 6
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    cp -puv --no-preserve=owner "${mummerout5a}1-0.png" "${mummerout5a}1-0.tn.jpg" "$nboutdir"/
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi

  if [ ! -s "${mummerout5b}1-0.tn.jpg" ]; then
    echo "5b. Comparing v3d $chr to v3c genome with MUMmer  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$v3cseq" \
      --outfile="$mummerout5b" \
      --size=large \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --savegnuplot \
      --nopoints \
      --minlength=5000 \
      --minpercent=96 \
      --showgaps=100000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    tn "${mummerout5b}1-0.png" 6
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    cp -puv --no-preserve=owner "${mummerout5b}1-0.png" "${mummerout5b}1-0.tn.jpg" "$nboutdir"/
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done  # for chr



# 6. Comparison of MUMmer organellar to self
for chr in MT PT ; do
  v3dchr="/vcru_share_s/carrot/dcarv3d/data/DCARv3d_${chr}.fna.gz"
  mummerout6a="$outdir/$pfx.mummer6a.$chr"
  mummerout6b="$outdir/$pfx.mummer6b.$chr"

  if [ ! -s "${mummerout6b}1-0.alignment.txt" ]; then
    echo "6a. Comparing $chr with MUMmer to self  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$v3dchr" \
      --outfile="$mummerout6a" \
      --size=large \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --savegnuplot \
      --nopoints \
      --showgaps=1000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    tn "${mummerout6a}1-0.png" 6
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
    cp -puv --no-preserve=owner "${mummerout6a}1-0.png" "${mummerout6a}1-0.tn.jpg" "$nboutdir"/
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi

    echo "6b. Comparing $chr with MUMmer to self with alignments  $(date)"
    bb.mummer --task=compare \
      --infile="$v3dchr" \
      --infile="$v3dchr" \
      --outfile="$mummerout6b" \
      --size=none \
      --savedeltafile --reusedeltafile \
      --savecoords \
      --showaligns \
      --nopoints \
      --showgaps=1000
    r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
  fi
done  # for chr



cp -puv --no-preserve=owner "$0" notebook/
echo "$0 Done  $(date)"
exit 0



#eof
