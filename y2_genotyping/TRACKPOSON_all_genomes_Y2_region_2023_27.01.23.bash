file=$1
file1=$(echo $file | awk -F "/" '{split($8,a,".");print a[1]}') 
fq1=$file
fq2=$( echo $file | sed -e s/_1/_2/)
out=$file1
#TE name
te="Y2_Helitron"
#path to genome_bed file
ref_genome_10kbpwindows="/home/Archive01/amacko/results_amp/TRACKPOSON/data/y2region_features_10kb_windows.bed"
#path to genome fasta file from NCBI
blast_ref_database="/home/Archive01/amacko/results_amp/TRACKPOSON/data/y2_region.fas"
#path to genome fasta file masked or no ins
blast_ref_database_no_ins="/home/Archive01/amacko/results_amp/TRACKPOSON/data/y2_region_no_hel.fas"
#bowtie2 index for TE
DB="/home/Archive01/amacko/results_amp/TRACKPOSON/data/Y2_Helitron"
#bowtie2 index for genome
DBG="/home/Archive01/amacko/results_amp/TRACKPOSON/data/y2_region"
#bowtie2 index for genome without insertion
DBGnoins="/home/Archive01/amacko/results_amp/TRACKPOSON/data/y2_region_no_hel"
#path to results for TE
DIR="/home/Archive01/amacko/results_amp/TRACKPOSON/results/Y2_region_2023"

cd $DIR

###### build bowtie2 index with the TE reference sequence
#bowtie2-build TE.fa TE


###### mapping reads against TE 
bowtie2 --time --end-to-end  -k 1 --very-fast -p 7 -x $DB -1 $fq1 -2 $fq2   | samtools view -bS -@ 7 - > "$out"-vs-"$te".bam


####### keep only unmap reads with flag unmap/map
# 69 = read paired,  read unmapped,  first in pair
# 133 = read paired,  read unmapped, second in pair
# 165 = read paired,  read unmapped, second in pair,mate reverse strand
# 181 =  read paired, read unmapped,  read reverse strand, mate reverse strand,second in pair
# 101 = read paired, read unmapped, mate reverse strand, first in pair
# 117 = read paired, read unmapped,   read reverse strand, mate reverse strand,   first in pair

	samtools view "$out"-vs-"$te".bam | awk -F "\t" '{if ( ($1!~/^@/) && (($2==69) || ($2==133) || ($2==165) || ($2==181) || ($2==101) || ($2==117)) ) {print ">"$1"\n"$10}}' > $out-vs-$te.fa


###### blast .fa against a reference genome for identification insertion
blastn -db $blast_ref_database -query $out-vs-$te.fa -out $out-vs-$te.fa.bl -num_threads 40 -evalue 1e-20

###### parse blast to find insertions

perl /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/find_insertion_point.pl $out-vs-$te.fa.bl $out-vs-$te

###### sort bed output
	sort -k1,1 -k2,2n $out-vs-$te.bed > $out-vs-$te.sort.bed
#

#####################################################################################################################
##### to identify soft-clipped reads; added to the original TRACKPOSON pipeline 
#####################################################################################################################
######  map against the genome masked/without TE - use all fq for bowtie2 mapping
bowtie2 --time  --very-fast-local  -p 40 -x $DBGnoins -1 $fq1 -2 $fq2  | samtools view -bS -@ 40 - > "$out"-vs-"$te".G_nohel_soft.bam

######## keep/print soft-clipped part of reads => part of TE
	/home/amacko/programy/SE-MEI/extractSoftclipped  "$out"-vs-"$te".G_nohel_soft.bam >"$out"-vs-"$te".Gsoft_nohel.fq.gz
	seqtk seq -a "$out"-vs-"$te".Gsoft_nohel.fq.gz > "$out"-vs-"$te".softnohel.fasta

#-------------------------------------------------------------------------------------------------------------------------------------
######  map against the genome with TE - use all fq for bowtie2 mapping
	bowtie2 --time  --very-fast-local  -p 40 -x $DBG -1 $fq1 -2 $fq2  | samtools view -bS -@ 40 - > "$out"-vs-"$te".Gsoft.bam

########keep/printsoft-clipped part of reads => part of TE
	/home/amacko/programy/SE-MEI/extractSoftclipped  "$out"-vs-"$te".Gsoft.bam >"$out"-vs-"$te".Gsoft.fq.gz
	seqtk seq -a "$out"-vs-"$te".Gsoft.fq.gz > "$out"-vs-"$te".softG.fasta

#-------------------------------------------------------------------------------------------------------------------------------------

#### blast split-reads
	
	blastn -db $blast_ref_database -query $out-vs-$te.softnohel.fasta -out $out-vs-$te.softnohel.bl -num_threads 40 -dust no -soft_masking false -word_size 10
	perl /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/find_insertion_point_soft.pl $out-vs-$te.softnohel.bl  $out-vs-$te.softnohel
	bedtools sort -i $out-vs-$te.softnohel.bed >$out-vs-$te.softnohel.sorted.bed 

	blastn -db $blast_ref_database -query  $out-vs-$te.softG.fasta -out $out-vs-$te.softG.bl -num_threads 40 -dust no -soft_masking false -word_size 10
	perl /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/find_insertion_point_soft.pl $out-vs-$te.softG.bl  $out-vs-$te.softG
	bedtools sort -i $out-vs-$te.softG.bed >$out-vs-$te.softG.sorted.bed 


#bedtools merge -i $out-vs-$te.soft.sorted.bed -c 1 -o count >$out-vs-$te.soft.merged.bed
######coveragebed by 10kb windows - to have higher cov - use 5 o higher number at this step
#bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.sort.bed | awk -F "\t" '{if ($4>=5){print $0}}' > coveragebed_5_$out-vs-$te\_per10kb.bed



##################### cov 

######coveragebed by 10kb windows - to have higher cov - use 5 o higher number at this step
	bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.sort.bed | awk -F "\t" '{if ($4>=1){print $0}}' > coveragebed_1_$out-vs-$te\_per10kb.bed

#10
	bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.sort.bed | awk -F "\t" '{if ($4>=3){print $0}}' > coveragebed_3_$out-vs-$te\_per10kb.bed


##################### cov soft

######coveragebed by 10kb windows - to have higher cov - use 5 o higher number at this step
	bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.softnohel.sorted.bed | awk -F "\t" '{if ($4>=1){print $0}}' > coveragebed_1_$out-vs-$te\_per10kb.soft.G_nohel.bed
#cov 3
	bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.softnohel.sorted.bed | awk -F "\t" '{if ($4>=3){print $0}}' > coveragebed_3_$out-vs-$te\_per10kb.soft.G_nohel.bed

######coveragebed by 10kb windows - to have higher cov - use 5 o higher number at this step
	bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.softG.sorted.bed | awk -F "\t" '{if ($4>=1){print $0}}' > coveragebed_1_$out-vs-$te\_per10kb.soft.G.bed
#cov 3
	bedtools coverage -counts -nonamecheck -a $ref_genome_10kbpwindows -b $out-vs-$te.softG.sorted.bed | awk -F "\t" '{if ($4>=3){print $0}}' > coveragebed_3_$out-vs-$te\_per10kb.soft.G.bed


######cleaning temporary files
rm $out-vs-$te.bam
rm $out-vs-$te.fa*
rm $out-vs-$te.bed





