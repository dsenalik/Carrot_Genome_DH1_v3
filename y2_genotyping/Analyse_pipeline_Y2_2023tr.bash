#####
#Date : 28.01.2023
# Y2_helitron
#####

#!/usr/bin/bash

#threshold 2 instead of threashold 5
# for file in *sort.bed;do out=$(echo $file | sed -e "s/\.sort/\_per10kb/");bedtools coverage -counts -b $file -a  ~/Bureau/Database/IRGSP-1.0_10kbpwindows.bed | awk  -F "\t" '{if ($4>=2){print $0}}' > coveragebed\_$out;done


#define TE name
te=Y2_Helitron
echo $te
mkdir final_cov1
cp coveragebed_1*per10kb.bed final_cov1
cd final_cov1

#get list of all ins sites
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#get in sites with threshold=2
for file in coveragebed_*;do awk -F "\t" '{if ($4>=2){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov1_$te.names
wc -l all_position_cov1_$te.names

#reformat data
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#clean evv/remove empty files
find ./ -size 0 -exec rm -f {} \;

#run script to get R matrix
R CMD BATCH "--args $te" /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/Analyse_pipeline_cov1.R

##reformat matrix

sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' matrice_final_$te.csv>matrice_final_"$te"_transposed.csv




cd ..

####### soft nohel - to get info supporting TE insertion ######
mkdir final_cov1snohel
cp coveragebed_1*per10kb.soft.G_nohel.bed final_cov1snohel
cd final_cov1snohel

#get list of all ins sites
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#get in sites with threshold=2
for file in coveragebed_*;do awk -F "\t" '{if ($4>=2){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov1_$te.names
wc -l all_position_cov1_$te.names

#reformat data
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#clean evv/remove empty files
find ./ -size 0 -exec rm -f {} \;

#run script to get R matrix
R CMD BATCH "--args $te" /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/Analyse_pipeline_cov1.R

###reformat matrix

sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' matrice_final_$te.csv>matrice_final_"$te"_transposed.csv

cd ..


####### soft nohel-to get info supporting TE insertion ######
mkdir final_cov1s
cp coveragebed_1*per10kb.soft.G.bed final_cov1s
cd final_cov1s

#get list of all ins sites
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#get in sites with threshold=2
for file in coveragebed_*;do awk -F "\t" '{if ($4>=2){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov1_$te.names
wc -l all_position_cov1_$te.names

#reformat data
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#clean evv/remove empty files
find ./ -size 0 -exec rm -f {} \;

#run script to get R matrix
R CMD BATCH "--args $te" /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/Analyse_pipeline_cov1.R

###reformat matrix

sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' matrice_final_$te.csv>matrice_final_"$te"_transposed.csv

cd ..





mkdir final_cov3
cp coveragebed_3*per10kb.bed final_cov3

#cp /a/transposon_M_CH_carrot/results/TAR_18/final_cov5/coveragebed_2*.bed final_cov2
cd final_cov3

#get list of all ins sites
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#get in sites with threshold=4
for file in coveragebed_*;do awk -F "\t" '{if ($4>=4){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov3_$te.names
wc -l all_position_cov3_$te.names

#reformat data
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#clean evv/remove empty files
find ./ -size 0 -exec rm -f {} \;

#run script to get R matrix
R CMD BATCH "--args $te" /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/Analyse_pipeline_cov3.R

###reformat matrix
sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' matrice_final_$te.csv>matrice_final_"$te"_transposed.csv

cd ..



mkdir final_cov3shel
cp coveragebed_3*per10kb.soft.G_nohel.bed final_cov3shel

#cp /a/transposon_M_CH_carrot/results/TAR_18/final_cov5/coveragebed_2*.bed final_cov2
cd final_cov3shel

#get list of all ins sites
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#get in sites with threshold=4
for file in coveragebed_*;do awk -F "\t" '{if ($4>=4){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov3_$te.names
wc -l all_position_cov3_$te.names

#reformat data
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#clean evv/remove empty files
find ./ -size 0 -exec rm -f {} \;

#run script to get R matrix
R CMD BATCH "--args $te" /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/Analyse_pipeline_cov3.R

###reformat matrix
sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' matrice_final_$te.csv>matrice_final_"$te"_transposed.csv

cd ..




mkdir final_cov3s
cp coveragebed_3*per10kb.soft.G.bed final_cov3s

#cp /a/transposon_M_CH_carrot/results/TAR_18/final_cov5/coveragebed_2*.bed final_cov2
cd final_cov3s

#get list of all ins sites
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#get in sites with threshold=4
for file in coveragebed_*;do awk -F "\t" '{if ($4>=4){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov3_$te.names
wc -l all_position_cov3_$te.names

#reformat data
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#clean evv/remove empty files
find ./ -size 0 -exec rm -f {} \;

#script R matrice
R CMD BATCH "--args $te" /home/Archive01/amacko/results_amp/TRACKPOSON/TRACKPOSON/Analyse_pipeline_cov3.R

###run script to get R matrix
sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' matrice_final_$te.csv>matrice_final_"$te"_transposed.csv

cd ..


