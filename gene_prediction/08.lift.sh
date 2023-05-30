#!/bin/bash
set -o pipefail
set -u # error if reference uninitialized variable
shopt -s nullglob

# lift annotations from v2 to v3b genome assemblies

# configuration
pfx=$(echo $0 | sed -e "s/^.*\///" -e "s/\..*$//")
chrs="PT"



# input files
ann="/gen/0061/2015cgp/scripts/15.DCARv2.sequintable"
inseqv2="/vcru_share_s/carrot/LNRQ01/data/130.DHv2_~.fna.gz"  # substitute ~ with chromosome
inseqv3="/vcru_share_s/carrot/dcarv3c/data/DCARv3c_~.fna.gz"



# output files
outdir="$pfx.outputdir"


# Initialization
if [ ! -d "$outdir" ]; then
  echo "Creating directory \"$outdir\"  $(date)"
  mkdir "$outdir"
  r=$?; if [ $r != 0 ]; then echo "Error code $r line $LINENO script $0, halting"; exit $r; fi
fi
cd "$outdir"



for chr in PT MT ; do
  echo "+++ Processing chromosome \"$chr\"  $(date)"
  v2=$(echo $inseqv2 | sed -e "s/~/$chr/")
  if [ ! -s "$v2" ]; then echo "Missing input file \"$v2\"" ; exit 1; fi
  v3=$(echo $inseqv3 | sed -e "s/~/$chr/")
  if [ ! -s "$v3" ]; then echo "Missing input file \"$v3\"" ; exit 1; fi
  sequin="$pfx.$chr.v2.sequin"
  inbed="$pfx.$chr.v2.bed"
  outbed="$pfx.$chr.v3.bed"
  outsequin="$pfx.$chr.v3.sequin"
  if [ ! -s "$sequin" ]; then
    echo "+++ Creating V2 sequin file for $chr"
    bb.fastagrep --infile="$ann" --outfile="$sequin" --grep="Feature DCARv2_$chr"
  fi

#set -x
set -e
  dir1=1.$chr.split
  dir2=2.$chr.lift
  dir3=3.$chr.psl
  dir4=4.$chr.liftup
  dir5=5.$chr.chain_raw
  dir6=6.$chr.chain_split
  dir7=7.$chr.net
  dir8=8.$chr.over
  
  mkdir -p $dir1
  mkdir -p $dir2
  mkdir -p $dir3
  mkdir -p $dir4
  mkdir -p $dir5
  # 6 made by chainSplit
  mkdir -p $dir7
  mkdir -p $dir8

  # uncompress fasta files temporarily
  echo "+++ Step 1 $chr Uncompress fasta files"
  zcat "$v2" > $pfx.$chr.old.tmp
  zcat "$v3" > $pfx.$chr.new.tmp
ls -ltr $pfx.$chr.*.tmp

  # split the new assembly sequences into 3K chunks and make lift files
  echo "+++ Step 2 $chr Split new assembly"
  faSplit size $pfx.$chr.new.tmp 3000 ./$dir1/$chr.split -lift=./$dir2/$chr.lft -oneFile
ls -ltr $dir1/
ls -ltr $dir2/
 
  # run blat, subject is old assembly, query is new assembly
  echo "+++ Step 3 $chr BLAT"
  blat "$pfx.$chr.old.tmp" ./$dir1/$chr.split.fa -t=dna -q=dna -tileSize=12 \
       -fastMap -minIdentity=95 -noHead -minScore=100 ./$dir3/$chr.psl
ls -ltr $dir3/

  # for plastid, we can ignore reverse matches due to the inverted repeat
  if [ "$chr" == "PT" ]; then
    mv ./$dir3/$chr.psl ./$dir3/$chr.psl.full
    grep -v '	-	' ./$dir3/$chr.psl.full > ./$dir3/$chr.psl
  fi

  # Change coordinates of .psl files to parent coordinate system
  echo "+++ Step 4 $chr Change coordinates"
  liftUp -pslQ ./$dir4/$chr.liftup.psl ./$dir2/$chr.lft warn ./$dir3/$chr.psl
ls -ltr $dir4/
 
  # Make chain files
  echo "+++ Step 5 $chr Made chain files"
  axtChain -linearGap=medium -faQ -faT -psl ./$dir4/$chr.liftup.psl \
           ./$pfx.$chr.old.tmp ./$pfx.$chr.new.tmp ./$dir5/$chr.chain
ls -ltr $dir5/
 
  # Merge and sort chain files
  echo "+++ Step 6a $chr Merged chain files"
  chainMergeSort ./$dir5/*.chain | chainSplit $dir6 stdin
ls -ltr $dir6/
 
  echo "+++ Step 6b $chr FASTA lengths"
  faSize $pfx.$chr.old.tmp -detailed > old.$chr.chr_length.txt
  faSize $pfx.$chr.new.tmp -detailed > new.$chr.chr_length.txt
ls -ltr *chr_length.txt
 
  # Make alignment nets from chain files
  for i in ./$dir6/*.chain
  do
  echo "+++ Step 7 $chr Process i=\"$i\""
    tag=${i/\.\/$dir6\//}
    chainNet $i ./old.$chr.chr_length.txt ./new.$chr.chr_length.txt ./$dir7/$tag.net /dev/null
  done
 
  # Create liftOver chain file
  for i in ./$dir6/*.chain
  do
  echo "+++ Step 8 $chr Process i=\"$i\""
    tag=${i/\.\/$dir6\//}
    netChainSubset ./$dir7/$tag.net $i ./$dir8/$tag.chain
  done
 
  echo "+++ Step 9 $chr chain file"
  cat ./$dir8/*.chain > $chr.old_new.over.chain
  ls -ltr $chr.old_new.over.chain

  # convert sequin to bed
  /gen/0081/scripts/06.sequintobed.pl "$sequin" "$inbed"
 
  # Make bed file to report converted coordinates. We can give the coordinates of
  # our query regions (based on old assembly) in the $inbed file and liftOver
  #  will report the converted coordinates in $outbed file.
  liftOver "$inbed" ./$chr.old_new.over.chain "$outbed" $chr.unMapped

  # convert bed back to sequin
  /gen/0081/scripts/06.sequintobed.pl "$outbed" "$outsequin"

  # remove temporary uncompressed fasta files
  rm -f $pfx.$chr.{old,new}.tmp

  echo "+++ Finished chromosome \"$chr\"  $(date)"

done  # for chr



cd -
cp -puv --no-preserve=owner "$0" notebook/
echo "$0 Done  $(date)"
exit 0



#eof
