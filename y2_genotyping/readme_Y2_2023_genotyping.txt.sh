
# Y2_region: DCARv3_Chr7:39164192..39186013 (+ strand) + 1kb (39163192-39187013)
# Y2_region_no_hel: DCARv3_Chr7:39164192..39186013 (+ strand) + 1kb (39163192-39187013) but helitron deleted from the sequence -> wild type


#1. Set env
######## configure env

#for thr soft-clipped reads:
 git clone https://github.com/dpryan79/SE-MEI.git
 cd SE-MEI
 git submodule update --init
 make

conda create --name TRACKPOSON -c daler sratoolkit
conda activate TRACKPOSON
conda install -c bioconda bowtie2 
conda install -c bioconda samtools 
conda install -c biocore blast-plus
conda install -c bioconda bedtools
conda install -c bioconda perl-bioperl
conda install -c bioconda picard
conda install -c bioconda htslib
conda deactivate TRACKPOSON



#2. Run modified TRACKPOSON 
conda activate TRACKPOSON


# Prepare genomic regions data:
#path/TRACKPOSON/data/Y2_Helitron.fas #TE
#path/TRACKPOSON/data/y2_region.fas # for transckposon
#path/TRACKPOSON/data/y2_region_no_hel.fas # for covarage at junction


cd path/TRACKPOSON/data/

bowtie2-build  Y2_Helitron.fas  Y2_Helitron #TE-fasta file in the folder with the index
bowtie2-build path/TRACKPOSON/data/y2_region.fas /home/Archive01/amacko/results_amp/TRACKPOSON/data/y2_region
bowtie2-build path/TRACKPOSON/data/y2_region_no_hel.fas /home/Archive01/amacko/results_amp/TRACKPOSON/data/y2_region_no_hel

#  blast+ database from the reference genome
makeblastdb -in path/TRACKPOSON/data/y2_region.fas -dbtype nucl -title y2_region
makeblastdb -in path/TRACKPOSON/data/y2_region_no_hel.fas -dbtype nucl -title y2_region_no_hel

#bwa index - genome index done
bwa index -p path/TRACKPOSON/data/bwa_y2 path/TRACKPOSON/data/y2_region.fas 
bwa index -p path/TRACKPOSON/data/bwa_y2nh path/TRACKPOSON/data/y2_region_no_hel.fas 
bwa index -p path/TRACKPOSON/data/Y2_Helitron path/TRACKPOSON/data/Y2_Helitron.fas


# Get a 10kb windows bed file from the reference genome 
samtools faidx y2_region.fas 
cut -f1,2 y2_region.fas.fai > size.y2region
bedtools makewindows  -g size.y2region -w 10000 > y2region_10kb_windows.bed # 3 bins 

# bed with bins representing helitron 3' end flanking regiin, helitron insertion, and helitron 3' flanking region and save as: y2region_features_10kb_windows.bed
#create a new dir for each varinat 

mkdir path/TRACKPOSON/results/Y2_region_2023 #for identification of insertion site
mkdir path/TRACKPOSON/results/Y2_region_no_hel_2023 #for identification of empty site
#than:TRACKPOSON/Analyse_pipeline.sh # 


##########################################
cd path/TRACKPOSON

for file in path/*1.fq.gz; do pathTRACKPOSON/TRACKPOSON_all_genomes_Y2_region_2023_27.01.23.bash $file;done

#1 1575
#2 3115 =>+1545 files (missing 35) - no reads, no insertion
#1545 extracted fasta
#blast 1575 - some empty
#1415 bed

#### 3. Run modified Analyse pipeline for the 

cd results/Y2_region_2023

tr -d '\r' <path/TRACKPOSON/Analyse_pipeline_Y2_2023.bash >path/TRACKPOSON/Analyse_pipeline_Y2_2023tr.bash

chmod +x  path/TRACKPOSON/Analyse_pipeline_Y2_2023tr.bash

path/TRACKPOSON/Analyse_pipeline_Y2_2023tr.bash 

