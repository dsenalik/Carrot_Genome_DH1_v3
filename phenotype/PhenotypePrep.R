# Step One: Prepare R Environment

# (1.0) Clear Global R 
rm(list = ls())
# (1.1) Set working directory ------------------------------> # Specific to computer
# (1.2) Load packages required for analysis
Pckg.Lst <-c("vroom","readxl","tidyverse", "PerformanceAnalytics",
             "lme4","lmerTest", "ggpubr", "corrplot")
package.check <- lapply(
  Pckg.Lst,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)}})
# (1.3) Load DataSheet from Ellison   
Raw.data <- read_csv("Raw.data.csv", col_names = T)  
View(Raw.data)

# Step Two: Identify outliers from the raw distribution of data.
# (2.1) Create Plot for Lutein
Raw.Lutein <- ggplot(data=Raw.data, aes(Raw.data$`ug/g dry Lutein`)) + 
  geom_histogram() + # Specify a ggplot - histogram
  xlab("Lutein (μg/g)") + # Update x axis label
  ylab("Observations") + # Change y axis label
  theme_classic() # Make the plot simple 
# (2.2) Create a histogram of the estimated concentration of each carotenoid
Raw.Alpha <- ggplot(data=Raw.data, aes(Raw.data$`ug/g dry Alpha`)) + 
  geom_histogram() +
  xlab("Alpha (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.3) Create a histogram of the estimated concentration of each carotenoid
Raw.Beta <- ggplot(data=Raw.data, aes(Raw.data$`ug/g dry Beta`)) + 
  geom_histogram() +
  xlab("Beta (μg/g)") +
  ylab("Observations") +
  theme_classic()
# (2.4) Smoosh Plots together for quick viewing
Raw.data.histograms <- ggarrange(Raw.Alpha,Raw.Beta,Raw.Lutein, 
          
                                 labels = c("A", "B", "C", "D"),
                                 ncol=2, nrow=2)
# (2.6) Clean Up the environment
rm(Raw.Alpha, Raw.Beta, Raw.Lutein) 
# (2.7) View Histogram  
Raw.data.histograms  
###  !!! There a samples with contrations that a much higher than other samples for Alpha, Beta, Lutein 
  ### !!!Fortunately, there are technical replicates for each sample where the technical replicates are the same sample injected in the HPLC twice.  
  #### The technical replicates can be used to find if both technical replicates are much higher.
Beta.High <- which(x=Raw.data$`ug/g dry Beta` > 2500)  
Raw.data$TechRep[Beta.High]
# (2.9) Identify samples with high Lutein
Lut.High <- which(x=Raw.data$`ug/g dry Lutein` > 300)  
Raw.data$TechRep[Lut.High]
# (2.10) Identify samples with high Alpha
Alpha.High <- which(x=Raw.data$`ug/g dry Alpha` > 2000)  
Raw.data$TechRep[Alpha.High]
# (2.11) # Best to check the lutein samples where both tech reps do not exceed threshold of 300
which(Raw.data$Sample.Name == "80209")
#[1] 325 326
Raw.data$`ug/g dry Lutein`[325:326]
# No reason to remove sample 80209 and tech replicate
which(Raw.data$Sample.Name == "80536")
# [1] 839 840
Raw.data$`ug/g dry Lutein`[839:840]
# Big difference... but will check the other samples before removing!

# Step Three: Identify any large batch effects in the HPLC data. 
  # Due to the number of samples not all HPLC data was collected at the same time. 
  # A standard was included to correct for "Batch" effects, in addition to the technical replicates. 
  # Technical replicates were completed on the same HPLC run with all the first technical replicates going first, followed by all of the second technical replicates. 
  # If there were issues with the HPLC-machine during the first or second technical replicate, they could be identifiable by consistently lower or higher estimated concentrations. 
  # To start identifying differences, I will look at the minimum, maximum and average of each HPLC batch concentration, and look differences in standard deviation between replicates. 
# (3.0) Identify each batch of HPLC completed by Standard
# Each run has a correction based on the Standard 
# This can be used to identify unique HPLC runs. 
Batches <- unique(Raw.data$Standard)
# Fifteen "batches' of HPLC completed! 
# (3.1) Create New Column to differentiate batches with "categorical variable"
Seasoned.Raw.Data <- Raw.data %>%
  mutate(Batches=Raw.data$Standard)
# (3.2) Create list replace Numeric Raw.Data$Standard 
Seasoning <- c("En", "To", "Tre", "Fire", "Fem",
               "Seks", "Sju", "Åtte", "Ni", "Ti",
               "Elleve", "Tolv", "Tretten", "Fjorten", "Femten")
# (3.3) Make new column to store categorical batch data
Seasoned.Raw.Data$Batches <- as.character(Seasoned.Raw.Data$Batches)
# (3.4) Replace Numeric batch information with categories  
for (i in 1:length(Batches)){
  Chef <- which(Seasoned.Raw.Data$Batches == Batches[i])
  Seasoned.Raw.Data$Batches[Chef] <- Seasoning[i]}

# (3.5) Check for differences in average, minimum, and maximum between batches. 
Batch.Check <- Seasoned.Raw.Data %>%
  group_by(Batches) %>%
  summarize(Alpha.Avg = mean(`ug/g dry Alpha`, na.rm=TRUE), 
            Alpha.Min = min(`ug/g dry Alpha`, na.rm=TRUE),
            Alpha.Max = max(`ug/g dry Alpha`, na.rm=TRUE),
            Beta.Avg = mean(`ug/g dry Beta`, na.rm=TRUE), 
            Beta.Min = min(`ug/g dry Beta`, na.rm=TRUE),
            Beta.Max = max(`ug/g dry Beta`, na.rm=TRUE),
            Lut.Avg = mean(`ug/g dry Lutein`, na.rm=TRUE), 
            Lut.Min = min(`ug/g dry Lutein`, na.rm=TRUE),
            Lut.Max = max(`ug/g dry Lutein`, na.rm=TRUE)) 
# (3.6) Look at results
View(Batch.Check)
# Hard to say off of this table. Plot the differences
# (3.7) Check distribution of samples!  
Alpha.Batch <- ggplot(data=Seasoned.Raw.Data,aes(x=Batches, y=`ug/g dry Alpha`, fill=Batches)) +
  geom_boxplot() + # Define boxplot
  scale_fill_manual(values=
                      rep(c("#C4012F", "White", "Black" ),5)) + # get UW colors
  theme_classic() + # keep it simple
  theme(legend.position = "none") # no legend needed
Beta.Batch <- ggplot(data=Seasoned.Raw.Data,aes(x=Batches, y=`ug/g dry Beta`, fill=Batches)) +
  geom_boxplot() +
  scale_fill_manual(values=
                      rep(c("#C4012F", "White", "Black" ),5)) +
  theme_classic() + 
  theme(legend.position = "none")
Lut.Batch <- ggplot(data=Seasoned.Raw.Data,aes(x=Batches, y=`ug/g dry Lutein`, fill=Batches)) +
  geom_boxplot() +
  scale_fill_manual(values=
                      rep(c("#C4012F", "White", "Black" ),5)) +
  theme_classic() + 
  theme(legend.position = "none")
# (3.8) Smoosh Plots together for quick viewing
Batch.Check.Boxes <- ggarrange(Alpha.Batch,Beta.Batch,Lut.Batch, 
                               ncol=1, nrow=4)
# (3.9) Clean Up the environment
rm(Alpha.Batch,Beta.Batch,Lut.Batch) 
# (3.10) View Histogram  
Batch.Check.Boxes  
#There are certainly differences in estimated concentration between HPLC runs. 
# However nothing is too egregious and as such no filtering is required after looking at these boxplots. 

# Step Four: Identify any inconsistency between technical replicates in the HPLC data.
# (4.1) Calculate Average Concentration 
# (4.1) Calculate Average Concentration 
Seared.Data <- Seasoned.Raw.Data %>%
  group_by(Sample.Name) %>%
  dplyr::summarize(`Alpha (ug/g)` = mean(`ug/g dry Alpha`, na.rm=TRUE), 
                   Alpha.sd = sd(`ug/g dry Alpha`, na.rm=TRUE),
                   `Beta (ug/g)` = mean(`ug/g dry Beta`, na.rm=TRUE), 
                   Beta.sd = sd(`ug/g dry Beta`, na.rm=TRUE),
                   `Lut (ug/g)` = mean(`ug/g dry Lutein`, na.rm=TRUE), 
                   Lut.sd = sd(`ug/g dry Lutein`, na.rm=TRUE))
# (4.2) Add Batch Information to compare the variation between technical replicates across batches. 
for (i in 1:length(Seared.Data$Sample.Name)){
  where <- match(Seared.Data$Sample.Name[i], Seasoned.Raw.Data$Sample.Name)
  Seared.Data[i,8] <- Seasoned.Raw.Data$Batches[where]
}
colnames(Seared.Data)[8] <- "Batch.Eff"
# (4.3) Plot Standard deviation between tech reps.  
Batch.Alpha <- ggplot(data=Seared.Data,aes(x=Batch.Eff, y=Alpha.sd)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  xlab("Standard Deviation Alpha") + # relabel X axix
  ylab("HPLC Batch") + # relabel Y axix
  theme_classic()  # Keep it simple
Batch.Beta <- ggplot(data=Seared.Data,aes(x=Batch.Eff, y=Beta.sd)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  xlab("Standard Deviation Beta") +
  ylab("HPLC Batch") +
  theme_classic()
Batch.Lutein <- ggplot(data=Seared.Data,aes(x=Batch.Eff, y=Lut.sd)) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  xlab("Standard Deviation Lutein") +
  ylab("HPLC Batch") +
  theme_classic()
# (4.4) Visualize differences
Var.Check <- ggarrange(Batch.Alpha, Batch.Beta, Batch.Lutein,
                       ncol=2, nrow=2)
# (4.5) Clean Up the environment
rm(Batch.Alpha, Batch.Beta, Batch.Lutein, ) 
# (4.6) View Histogram  
Var.Check  
# It appears there are a few samples with large variation between tech replicates for:
  # Tolv: Alpha, Beta
  # Ni Lutein
  # Elleve: Something wrong with this HPLC run. 

# Step Five: Find and remove Outliers in "Tolv (12th)" HPLC run
# (5.1) Subset a dataframe to focus on the 'Tolv" HPLC run
Tolv.dat <- subset(Seared.Data, Batch.Eff == "Tolv")
# (5.2) Identify outliers 
High.Alpha <- which(Tolv.dat$Alpha.sd > 200)
High.Beta <-which(Tolv.dat$Beta.sd > 800)
# (5.3) Identify outliers 
To.much.var <- c(High.Alpha, High.Beta)
# The 26th (Beta), 27th (All three), and 28th (Alpha, Phyto) were the problems
# (5.4) Identify outliers
To.much.var <- unique(To.much.var)
# (5.5) Get Sample Names 
Sour.Sample <- Tolv.dat$Sample.Name[To.much.var]
# (5.6) Remove HPLC data for the problem samples 
Out.Tech.Reps <- which(Seared.Data$Sample.Name == Sour.Sample)
Seared.Data[Out.Tech.Reps,2:8] <- NA
# (5.7) Clean R environment  
rm(High.Alpha,High.Beta,Tolv.dat, Sour.Sample)
# Samples 80036, 80037, 80038 had large differences between technical replicates. 
  # The fact the the samples before and after these three samples could indicate a pipetting issue rather than an HPLC problem. 
  # Regardless the HPLC data for these samples will be set to NA.

# Step Six: Find and remove Outliers in "Ni (9th)" HPLC run
# (6.1) Subset a dataframe to focus on the 'Ni" HPLC run
Ni.dat <- subset(Seared.Data, Batch.Eff == "Ni")
# (6.2) Identify outliers 
High.Lut <- which(Ni.dat$Alpha.sd > 50)
High.Lyco <-which(Ni.dat$Beta.sd > 50)
# (6.3) Identify outliers 
To.much.var <- c(High.Lut, High.Lyco)
# The 16th sample is the outlier
# (6.4) Identify outliers
To.much.var <- unique(To.much.var)
# (6.5) Get Sample Names 
Sour.Sample <- Ni.dat$Sample.Name[To.much.var]
# (6.6) Remove HPLC data for the problem samples 
Out.Tech.Reps <- which(Seared.Data$Sample.Name == Sour.Sample)
Seared.Data[Out.Tech.Reps,2:8] <-NA
# (6.7) Clean R environment  
rm(High.Lut, High.Lyco, Ni.dat, To.much.var, Sour.Sample)
# Sample 80432 had a large differences between technical replicates. 
# The fact the the samples before and after these three samples could indicate a pipetting issue rather than an HPLC problem. 
# The HPLC data for these samples will be set to NA.

# Step Seven: Figure out the problem with the Elleve (11th) HPLC Run
# (7.0) Subset to look at Elleve.dat
Elleve.dat <- subset(Seared.Data, Batch.Eff == "Elleve")
# (7.1) Check if there are consistent differences across tech replciates 
# Start by subsetting the Raw.data
where <- which(Seasoned.Raw.Data$Sample.Name %in% Elleve.dat$Sample.Name)
Seasoned.Elleve.dat <- Seasoned.Raw.Data[where,]

# (7.2) Sort the identify second tech replicates by 'B" added to second tech rep
her <- grepl("B", Seasoned.Elleve.dat$TechRep)
where <- which(her==TRUE)
# (7.3) Make matrix to store information about technical replicates
problem <- matrix(nrow=4, ncol=4)
problem[2:4,1] <- c("Avg.TRep1", "Avg.TRep2", "%Diff")  
problem[1,2:4] <- c("Alpha", "Beta", "Lut")
# (7.4) Calculate Averages
# Alpha
problem[2,2] <- mean(Seasoned.Elleve.dat$`ug/g dry Alpha`[-where])
problem[3,2] <- mean(Seasoned.Elleve.dat$`ug/g dry Alpha`[where])
problem[4,2] <- as.numeric(problem[3,2])/as.numeric(problem[2,2])
# Beta
problem[2,3] <- mean(Seasoned.Elleve.dat$`ug/g dry Beta`[-where])
problem[3,3] <- mean(Seasoned.Elleve.dat$`ug/g dry Beta`[where])
problem[4,3] <- as.numeric(problem[3,3])/as.numeric(problem[2,3])
# Lutein
problem[2,4] <- mean(Seasoned.Elleve.dat$`ug/g dry Lutein`[-where])
problem[3,4] <- mean(Seasoned.Elleve.dat$`ug/g dry Lutein`[where])
problem[4,4] <- as.numeric(problem[3,4])/as.numeric(problem[2,4])
# (7.5) Check out the results
View(problem)
# (7.6) Store Output
write.csv(x=problem, file="Summary.of.problem.HPLC.csv", row.names = F)
# There is a consistent decrease from technical replciate one to technical replcate 2. 
  # It is hard to determine when the HPLC results began to reduce. 
  # These samples set to NA 
# (7.7) Remove those samples in the Elleve HPLC run (Step 7) 
Problem.Samples <- which(Seared.Data$Sample.Name %in% Seasoned.Elleve.dat$Sample.Name)
Seared.Data[Problem.Samples,2:8] <- NA

# Step Eight: Color to HPLC reconciling 
  # Each sample has a categorical color score of Orange, Purple, Red, Yellow, and White. 
  # Yellow scores are expected to have high Lutein, 
  # Orange high Alpha and Beta. 
  # Knowing these we can check the HPLC vs color scores to make sure samples have appropriate values. 
  # The color scores were visually scored and entered into a spreadsheet. 
  #  It can be difficult to determine the difference between yellow and orange, 
  # as such some samples that were borderline were set to NA. 
  # Likewise there are samples that are orange, but have little to no Beta-Carotene in the HPLC data. 
  # The HPLC data for these samples will be set to NA. 
# (8.1) Load HPLC color discrepancy list
Color.HPLC.discrep <- read.table(file="HPLCQuestionableList.txt", header=F)
# (8.2) Remove Orange carrots where no Carotene was detected in HPLC
Problem.Samples.Color <- which(Seared.Data$Sample.Name %in% Color.HPLC.discrep$V1)
Seared.Data[Problem.Samples.Color,2:8] <- NA

# (8.3) Add color scores to seared data!
Medium.Rare.Data <-  Seared.Data
Colors <- read_csv("Color.Scores.csv", col_names=T)
colnames(Colors) <- c("Sample.Name", "Color.Score")
Colors <- Colors %>%
  filter(Sample.Name %in% Medium.Rare.Data$Sample.Name)
Medium.Rare.Data <- Medium.Rare.Data %>%
  filter(Sample.Name %in% Colors$Sample.Name)
Medium.Rare.Data <- Medium.Rare.Data[,-8]
Medium.Rare.Data$Sample.Name <- as.character(Medium.Rare.Data$Sample.Name)
Medium.Rare.Data <-left_join(Medium.Rare.Data,Colors, by="Sample.Name")
#(8.4) Convert Colors to a Numerica
Orange <- which(Medium.Rare.Data$Color.Score == "Orange")
Medium.Rare.Data$Color.Score[Orange] <- 1
Yellow <- which(Medium.Rare.Data$Color.Score == "Yellow")
Medium.Rare.Data$Color.Score[Yellow] <- 2
White <- which(Medium.Rare.Data$Color.Score == "White")
Medium.Rare.Data$Color.Score[White] <- 3
Red <- which(Medium.Rare.Data$Color.Score == "Red")
Medium.Rare.Data$Color.Score[Red] <- 4
Purple <- which(Medium.Rare.Data$Color.Score == "Purple")
Medium.Rare.Data$Color.Score[Purple] <- 5

######################################################################
# (9.0) Drop Columns that won't be used in GWA
Well.Done.Data <- select(Medium.Rare.Data, -Alpha.sd, -Beta.sd,
                         -Lut.sd)
# (9.1) Calculate Total Carotenoids
Well.Done.Data <- Well.Done.Data %>%
  mutate(`Total.Carotenoid(ug/g)` = `Alpha (ug/g)` + `Beta (ug/g)`+
           `Lut (ug/g)`)
# (9.2) Calculate Percentage of Carotenoids
Well.Done.Data <- Well.Done.Data %>%
  mutate(`Ratio Alpha` = `Alpha (ug/g)`/`Total.Carotenoid(ug/g)`)
Well.Done.Data <- Well.Done.Data %>%
  mutate(`Ratio Beta` = `Beta (ug/g)`/`Total.Carotenoid(ug/g)`)
Well.Done.Data <- Well.Done.Data %>%
  mutate(`Ratio Lutein` = `Lut (ug/g)`/`Total.Carotenoid(ug/g)`)
Well.Done.Data <- Well.Done.Data %>%
  mutate(`Ratio alpha_beta` = `Ratio Alpha` + `Ratio Beta`)
# (9.3) Write the phenotypic data to file
write.csv(x=Well.Done.Data, row.names = F, file="GWA.Phenotype.csv")


