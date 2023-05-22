# Population genomics identifies genetic signatures of carrot domestication and improvement and uncovers the origin of high carotenoid orange carrots

This repository contains scripts used for the following publication:
~Placeholder for reference to publication~



## Genome Assembly

### Genome Assembly Phase I

* FALCON (v.0.3.0) \
source: [https://github.com/PacificBiosciences/FALCON](https://github.com/PacificBiosciences/FALCON) \
parameters: `length_cutoff=5000 length_cutoff_pr=1000 --max_diff=50 --max_cov=50 --min_cov=4` \
reference: Chin, C.-S. et al. Phased diploid genome assembly with single-molecule real-time sequencing. Nature Methods 13, 1050-1054, [doi:10.1038/nmeth.4035](https://doi.org/10.1038/nmeth.4035) (2016).

* CANU (V.1.5) \
source: [https://github.com/marbl/canu/releases](https://github.com/marbl/canu/releases) \
reference: Koren, S. et al. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Res 27, 722-736, [doi:10.1101/gr.215087.116](https://doi.org/10.1101/gr.215087.116) (2017).

### Polishing

* CANU (V.1.5) \
source: [https://github.com/marbl/canu/releases](https://github.com/marbl/canu/releases) \
reference: Koren, S. et al. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Res 27, 722-736, [doi:10.1101/gr.215087.116](https://doi.org/10.1101/gr.215087.116) (2017).

* Arrow as implemented in GenomicConsensus (v?) \
source: [https://github.com/PacificBiosciences/GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)

* Pilon (V?) \
parameters: using default parameters \
reference: Walker, B. J. et al. Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLOS ONE 9, e112963, [doi:10.1371/journal.pone.0112963](https://doi.com/10.1371/journal.pone.0112963) (2014).

### Scaffolding

* BWA mem (0.7.17) \
reference: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754-1760, [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) (2009).

* SALSA2 (V?)
200 iterations and GATC as the restriction site \
reference: Ghurye, J. et al. Integrating Hi-C links with assembly graphs for chromosome-scale assembly. PLOS Computational Biology 15, e1007273, [doi:10.1371/journal.pcbi.1007273](https://doi.org/10.1371/journal.pcbi.1007273) (2019).

### Gap Fill

* PBJelly (V?) \
parameters: `--minMatch 8 --minPctIdentity 80 --bestn 1 --nCandidates 20 --maxScore -500 --nproc 39 --noSplitSubreads` \
reference: English, A. C. et al. Mind the gap: upgrading genomes with Pacific Biosciences RS long-read sequencing technology. PLoS One 7, e47768, [doi:10.1371/journal.pone.0047768](https://doi.org/10.1371/journal.pone.0047768) (2012).

### Genome Assembly Phase II

* BWA mem (0.7.17)
filtered for MAPQ>=10 and marker sequences with mapping frequency of >1 were discarded \
reference: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754-1760, [ [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) (2009).

* STAR (V?)
SAM files were then filtered for MAPQ >=10 and the reads that had both pairs mapping were kept \
reference: Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21, [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635) (2012).

* BWA mem (0.7.17) \
parameters: `-5 -t 44`

* sort \
parameters: `-k3,3d -k7,7d`

* juicer_tools.1.8.9_jcuda.0.8.jar \
parameters: `pre -q 10` \
reference: Durand, N. C. et al. Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments. Cell Syst 3, 95-98, [doi:10.1016/j.cels.2016.07.002](https://doi.org/10.1016/j.cels.2016.07.002) (2016).

* A custom program was used to visualize and identify connections of PE sequences (PE BACs, 10, 20 and 40 kb MPE) to other scaffolds or contigs


### Genome Assembly Phase III

* RaGOO (V?) \
parameters: `-C -t 43 -i 0.3` \
reference: Alonge, M. et al. RaGOO: fast and accurate reference-guided scaffolding of draft genomes. Genome Biology 20, 224, [doi:10.1186/s13059-019-1829-6](https://doi.org/10.1186/s13059-019-1829-6) (2019).

* LR-Gapcloser (V?) \
parameters: `-s p -t 40 -r 5` \
reference: Xu, G.-C. et al. LR_Gapcloser: a tiling path-based gap closer that uses long reads to complete genome assembly. GigaScience 8, [doi:10.1093/gigascience/giy157](https://doi.org/10.1093/gigascience/giy157) (2018).



## Assembly Quality Verification

### Sequence contamination assessment

* Fastq-Screen (V?) \
reference: Wingett, S. & Andrews, S. FastQ Screen: A tool for multi-genome mapping and quality control [version 2; peer review: 4 approved]. F1000Research 7, [doi:10.12688/f1000research.15931.2](https://doi.org/10.12688/f1000research.15931.2) (2018).

* GC content distribution estimates (?)

### Linkage map marker mapping

* BWA mem (0.7.17) \
reference: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754-1760, [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) (2009).

### Gene space coverage

* BWA mem (EST mapping) (0.7.17) \
source for EST sequences: [https://www.carrotomics.org/bio_data/291864](https://www.carrotomics.org/bio_data/291864)

* StringTie (transcriptome mapping) (V?) \
reference: Pertea, M. et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology 33, 290-295, [doi:10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122) (2015).

* GMAP (IsoSeq mapping) (V?) \
reference: Wu, T. D. & Watanabe, C. K. GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21, 1859-1875, [doi:10.1093/bioinformatics/bti310](https://doi.org/10.1093/bioinformatics/bti310) (2005).



## Repetitive Sequence Annotation

### De Novo detection

* RepeatModeler (v.2.0.1) \
source: [http://www.repeatmasker.org/RepeatModeler/](http://www.repeatmasker.org/RepeatModeler/)

### Masking

* RepeatMasker (v.4.1.0) \
source: [http://www.repeatmasker.org](http://www.repeatmasker.org)

### Identification, annotation, and age analysis of LTR-RTs

* LTRharvest (GenomeTools v.1.6.1) \
parameters: `-seed 80 -maxlenltr 4000 -mindistltr 3000 -mintsd 2 -maxtsd 20 -motif tgca` \
reference: Ellinghaus, D., Kurtz, S. & Willhoeft, U. LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC Bioinformatics 9, 18, [doi:10.1186/1471-2105-9-18](https://doi.org/10.1186/1471-2105-9-18) (2008).

* SiLiX (v1.2.9) \
parameters: `--ident 0.6 --overlap 0.7` \
reference: Miele, V., Penel, S. & Duret, L. Ultra-fast sequence clustering from similarity networks with SiLiX. BMC Bioinformatics 12, 116, [https://doi.org/10.1186/1471-2105-12-116](doi:10.1186/1471-2105-12-116) (2011).

* DANTE (v.1.1.0)

* Mafft (v.7.471) \
parameters: default settings \
reference: Katoh, K. & Standley, D. M. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol 30, 772-780, [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010) (2013).

* R ‘ape’ package (v.5.4-1) \
parameters: “K80” model \
reference: Paradis, E. & Schliep, K. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35, 526-528, [doi:10.1093/bioinformatics/bty633](https://doi.org/10.1093/bioinformatics/bty633) (2018).

* HelitronScanner \
parameters: default settings \
reference: Xiong, W., He, L., Lai, J., Dooner, H. K. & Du, C. HelitronScanner uncovers a large overlooked cache of <i>Helitron</i> transposons in many plant genomes. Proceedings of the National Academy of Sciences 111, 10263-10268, [doi:10.1073/pnas.1410068111](https://doi.org/10.1073/pnas.1410068111) (2014).

### Quality of the assembled repetitive sequences

* LTR Assembly Index (LAI) \
reference: Ou, S., Chen, J. & Jiang, N. Assessing genome assembly quality using the LTR Assembly Index (LAI). Nucleic Acids Research 46, e126-e126, [doi:10.1093/nar/gky730](https://doi.org/10.1093/nar/gky730) (2018).



## Gene Prediction and Genome Annotation

### RNA-seq from 20 DH1 tissues

* STAR (V?) \
reference: Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21, [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635) (2012).

* StringTie (V?) \
source: [https://ccb.jhu.edu/software/stringtie/](https://ccb.jhu.edu/software/stringtie/) \
reference: Pertea, M. et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology 33, 290-295, [doi:10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122) (2015).

### PacBio IsoSeq transcripts

* IsoSeq3 (v?) \
source: [https://github.com/PacificBiosciences/IsoSeq](https://github.com/PacificBiosciences/IsoSeq)

* GMAP (v?) \
parameters: -f=2; and minimum identity and coverage=99% \
reference: Wu, T. D. & Watanabe, C. K. GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21, 1859-1875, [doi:10.1093/bioinformatics/bti310](https://doi.org/10.1093/bioinformatics/bti310) (2005).

### Gene Model Prediction with MAKER

* MAKER (v?) \
reference: Holt, C. & Yandell, M. MAKER2: an annotation pipeline and genome-database management tool for second-generation genome projects. BMC Bioinformatics 12, 491, [doi:10.1186/1471-2105-12-491](https://doi.org/10.1186/1471-2105-12-491) (2011).

* AUGUSTUS (v2.5.5)

* SNAP (v?) \
source: [https://github.com/KorfLab/SNAP](https://github.com/KorfLab/SNAP)

### Gene Model Prediction with GeMoMa

* GeMoMa (v?) \
reference: Keilwagen, J., Hartung, F. & Grau, J. GeMoMa: Homology-Based Gene Prediction Utilizing Intron Position Conservation and RNA-seq Data. Methods Mol Biol 1962, 161-177, [doi:10.1007/978-1-4939-9173-0_9](https://doi.org/10.1007/978-1-4939-9173-0_9) (2019).

### Gene Model Curation

* GMAP (v?) \
reference: Wu, T. D. & Watanabe, C. K. GMAP: a genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics 21, 1859-1875, [doi:10.1093/bioinformatics/bti310](https://doi.org/10.1093/bioinformatics/bti310) (2005).

* GenomeThreader (v?) \
Reference: Gremme, G., Brendel, V., Sparks, M. E. & Kurtz, S. Engineering a software tool for gene structure prediction in higher organisms. Information and Software Technology 47, 965-978, [doi:https://doi.org/10.1016/j.infsof.2005.09.005](https://doi.org/10.1016/j.infsof.2005.09.005) (2005).

### Gene annotation and quality verification

* Blast2Go (v?) \
reference: Conesa, A. et al. Blast2GO: a universal tool for annotation, visualization and analysis in functional genomics research. Bioinformatics 21, 3674-3676, [doi:10.1093/bioinformatics/bti610](https://doi.org/10.1093/bioinformatics/bti610) (2005).

* PlantTFcat (v?) \
reference: Dai, X., Sinharoy, S., Udvardi, M. & Zhao, P. X. PlantTFcat: an online plant transcription factor and transcriptional regulator categorization and analysis tool. BMC Bioinformatics 14, 321, [doi:10.1186/1471-2105-14-321](https://doi.org/10.1186/1471-2105-14-321) (2013).

* PRGdb (v3.0) \
reference: Osuna-Cruz, C. M. et al. PRGdb 3.0: a comprehensive platform for prediction and analysis of plant disease resistance genes. Nucleic Acids Res 46, D1197-d1201, [doi:10.1093/nar/gkx1119](https://doi.org/10.1093/nar/gkx1119) (2018).

### Non-coding RNAs

* INFERNAL (v1.1.2) \
parameters: Rfam database (v13.0) and the p-value of <0.05 \
reference: Nawrocki, E. P. & Eddy, S. R. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29, 2933-2935, [doi:10.1093/bioinformatics/btt509](https://doi.org/10.1093/bioinformatics/btt509) (2013).

### Isoform Analysis

* Sqanti3 (v?) \
reference: Tardaguila, M. et al. SQANTI: extensive characterization of long-read transcript sequences for quality control in full-length transcriptome identification and quantification. Genome Res 28, 396-411, [doi:10.1101/gr.222976.117](https://doi.org/10.1101/gr.222976.117) (2018).



## Resequencing and Phenotyping

### Phenotyping

* Scripts and datafiles related to phenotyping can be found in [/phenotype/](/phenotype/)



## Variant Detection

### Mapping

* BWA (v 0.7.17-r1188) \
parameters: `-a -M -t 2 -R` readgroup \
reference: Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics 25, 1754-1760, [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) (2009).

* The Genome Analysis Toolkit (GATK) v4.0.7.0 \
example script for one genotype: [variant_detection/job1800.sh](variant_detection/job1800.sh) \
example joint genotyping script for one region: [variant_detection/template.jointgenodcarv3d.sh](variant_detection/template.jointgenodcarv3d.sh) \
reference: McKenna, A. et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20, 1297-1303, [doi:10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110) (2010).

### Genotyping the Y2 insertion

* A heavily modified version of [TRACKPOSON](https://github.com/Markestine/TRACKPOSON) \
scripts: [/y2_genotyping/](/y2_genotyping/)



## Population structure, phylogenetic, PCA and gene flow analysis

### Population Genetic Ancestry

* Plink (v1.90b3.44) \
parameters: `--indep 50 5 2` \
reference: Purcell, S. et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81, 559-575, [doi:10.1086/519795](https://doi.org/10.1086/519795) (2007).

* PHYLIP (v?) \
reference: Felsenstein, J. PHYLIP (phylogeny inference package), version 3.5 c.  (Joseph Felsenstein., 1993).

### Principal Component Analysis

* snpgdsPCA implemented in the R package SNPRelate (v1.20.1) \
Zheng, X. et al. A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics 28, 3326-3328, [doi:10.1093/bioinformatics/bts606](https://doi.org/10.1093/bioinformatics/bts606) (2012).

### Gene Flow

* TreeMix (v1.13) \
reference: Fitak, R. R. OptM: estimating the optimal number of migration edges on population trees using Treemix. Biology Methods and Protocols 6, bpab017, [doi:10.1093/biomethods/bpab017](https://doi.org/doi:10.1093/biomethods/bpab017) (2021).



## Genetic diversity, FST, linkage disequilibrium and demographic analysis

### Pairwise FST and nucleotide diversity (pi)

* VCFtools (v?) \
reference: Weir, B. S. & Cockerham, C. C. Estimating F-Statistics for the Analysis of Population Structure. Evolution 38, 1358-1370, [doi:10.2307/2408641](https://doi.org/10.2307/2408641) (1984).

### LD Decay

* PopLDdecay (v3.31) \
parameters: `-OutStat` \
reference: Zhang, C., Dong, S. S., Xu, J. Y., He, W. M. & Yang, T. L. PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. Bioinformatics 35, 1786-1788, [doi:10.1093/bioinformatics/bty875](https://doi.org/10.1093/bioinformatics/bty875) (2019).

### Estimates of effective population size history and divergence times

* bcftools (v1.10.2) \
parameters: `view -e 'F_MISSING >= 0.1' -Oz` \
reference: Danecek, P. et al. Twelve years of SAMtools and BCFtools. GigaScience 10, giab008, [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) (2021).

* SMC++ software (v.1.15.2) \
parameters and code: Code for the estimation of effective population size and divergence times using SMC++ is available at [https://github.com/mishaploid/carrot-demography](https://github.com/mishaploid/carrot-demography). \
source: [https://github.com/popgenmethods/smcpp](https://github.com/popgenmethods/smcpp) \
reference: Terhorst, J., Kamm, J. A. & Song, Y. S. Robust and scalable inference of population history from hundreds of unphased whole genomes. Nat Genet 49, 303-309, [doi:10.1038/ng.3748](https://doi.org/10.1038/ng.3748) (2017).



## Genome-Wide Scans for Signatures of Selection

### Region Identification

* VCFtools (v0.1.14 120) \
reference: Danecek, P. et al. The variant call format and VCFtools. Bioinformatics (Oxford, England) 27, 2156-2158, [doi:10.1093/bioinformatics/btr330](https://doi.org/10.1093/bioinformatics/btr330) (2011).

* XP-CLR (v1.0) \
reference: Chen, H., Patterson, N. & Reich, D. Population differentiation as a test for selective sweeps. Genome Research 20, 393-402 (2010).



## GWA analysis

### GWA Analysis
Scripts and parameters for this section are here: [/genotype/](/genotype/)

* vcftools (v0.1.16) \
reference: Weir, B. S. & Cockerham, C. C. Estimating F-Statistics for the Analysis of Population Structure. Evolution 38, 1358-1370, [doi:10.2307/2408641](https://doi.org/10.2307/2408641) (1984).

* Beagle (v5.0) \
parameters: default settings \
reference: Browning, S. R. & Browning, B. L. Rapid and Accurate Haplotype Phasing and Missing-Data Inference for Whole-Genome Association Studies By Use of Localized Haplotype Clustering. The American Journal of Human Genetics 81, 1084-1097, [doi:10.1086/521987](https://doi.org/10.1086/521987) (2007).

* Tassel (v5) \
reference: Bradbury, P. J. et al. TASSEL: software for association mapping of complex traits in diverse samples. Bioinformatics 23, 2633-2635, [doi:10.1093/bioinformatics/btm308](https://doi.org/10.1093/bioinformatics/btm308) (2007).

* GAPIT (custom modified version of 2020.10.24 Gapit 3.0) \
Script: [GWA/GWA.R.txt](GWA/GWA.R.txt) \
reference: Lipka, A. E. et al. GAPIT: genome association and prediction integrated tool. Bioinformatics 28, 2397-2399, [doi:10.1093/bioinformatics/bts444](https://doi.org/10.1093/bioinformatics/bts444) (2012).

### Significance Threshold

* simpleM (version July 2, 2009) \
script: [GWA/simpleMSigThreshold.R.txt](GWA/simpleMSigThreshold.R.txt) \
reference: Gao, X., Starmer, J. & Martin, E. R. A multiple testing correction method for genetic association studies using correlated single nucleotide polymorphisms. Genet Epidemiol 32, 361-369, [doi:10.1002/gepi.20310](https://doi.org/10.1002/gepi.20310) (2008).

### Manhattan Plots

* qqman (v0.1.8) \
script: [GWA/qqman4Manhattanplots.R](GWA/qqman4Manhattanplots.R) \
reference: Turner, S. qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots.  (2014).



## RNAseq analysis for Or and Y2

### transcriptome profile of candidate genes underlying the Or and Y2 locus

* TRIMMOMATIC (v0.36) \
reference: Bolger, A. M., Lohse, M. & Usadel, B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30, 2114-2120, [doi:10.1093/bioinformatics/btu170](https://doi.org/10.1093/bioinformatics/btu170) (2014).

* Rsubread (v?) in R 3.5.0 \
reference: Liao, Y., Smyth, G. K. & Shi, W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. Nucleic Acids Res 47, e47, [doi:10.1093/nar/gkz114](https://doi.org/10.1093/nar/gkz114) (2019).

* FeatureCounts (v?) in R 3.5.0 \
reference: Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923-930, [doi:10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656) (2014).



## Genetic effect and interaction analysis

### Alternative genetic effects including additive, dominance, recessive and over-dominance

* lm in R 3.5.0 \
reference: Everitt, B. Book reviews : Chambers JM, Hastie TJ eds 1992: Statisti cal models in S. California: Wadsworth and Brooks/Cole. ISBN 0 534 16765-9. Statistical Methods in Medical Research 1, 220-221, [doi:10.1177/096228029200100208](https://doi.org/10.1177/096228029200100208) (1992).
