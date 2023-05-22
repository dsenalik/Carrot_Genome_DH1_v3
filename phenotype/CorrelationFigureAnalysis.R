# Step One: Prepare R Environment
# (1.0) Load The TidyVerse 
  library('tidyverse')
# (2.0) Load Phenotypic data
  pheno <- read_csv("CorrelationFigureData.csv", col_names = T)
    colnames(pheno)[5] <- "Color"
    colnames(pheno)[3] <- "Carotene"
    colnames(pheno)[6] <- "Total.Carotenoid(ug/g)"
  pheno$Color <- as.factor(pheno$Color)
# (3.0) Compare alpha and beta in different color categories
  FigureB1b <-  ggplot(data=pheno, aes(x=Color, y=Carotene, fill=Color)) +
    geom_boxplot() + # Define boxplot
    scale_fill_manual(values= c("Dark Orange", "#FFA500", "Yellow", "White")) +
    ylab("α + β Carotene / Total Carotenoids (μg/g)") +
    theme_classic() + # keep it simple
    theme(legend.position = "none")
# (4.0) Test correlation between total carotenoids and ratio of Carotene
  improvement <- cor.test(pheno$Carotene, pheno$`Total.Carotenoid(ug/g)`, 
                          method="spearman")
# (5.0) Extract P value and R value from correlation test
  R <- improvement$estimate
  p <- improvement$p.value
# (5.1) plot the correlation
  FigureB1a <-  pheno %>%
    ggplot(aes(x=`Total.Carotenoid(ug/g)`, y=Carotene)) +
    geom_point() +
    geom_smooth(method="loess", mapping = aes(col="Red"), show.legend = F) +
    xlab("Total Carotenoids (μg/g)") +
    ylab("α + β Carotene / Total Carotenoids (μg/g)") +
    theme_classic()
# (5.2) Add Correlation r and p values to plot 
  FigureB1a <- FigureB1a  + annotate(geom="text", x=2000, y=0.3, label="r = 0.6677",
                                   color="black") + 
    annotate(geom="text", x=2000, y=0.2, label="p < 2.2e-16",
           color="black")
# (5.3) Combine  the figures
  library(ggpubr)
  png(filename="FigureB1.png", width=960, height=480, bg="white")
  ggarrange(FigureB1a, FigureB1b, ncol=2, nrow=1)
  dev.off()
  
  