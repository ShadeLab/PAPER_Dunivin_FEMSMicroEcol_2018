wd <- paste(getwd())
stats <- read_delim(paste(wd, "/../assembly_assessments/summary_stats_final.txt", sep = ""), delim = " ", col_types = list(col_character(),col_character(), col_number(),col_number(),col_number(),col_number(),col_number(),col_number(),col_number()))

#remove all extra columns
stats <- stats[!stats$ProteinContigClusters.99 == "ProteinContigClusters.99",]
stats$cen01_acr3_stats.txt <- gsub("arsC_", "arsC", stats$cen01_acr3_stats.txt)
stats$Site <- gsub("cen", "Cen", stats$cen01_acr3_stats.txt)

#tidy stats data
stats.tidy <- stats %>%
  separate(Site, into = c("Site", "Gene")) %>%
  left_join(meta, by = "Site")

ggplot(stats.tidy, aes(x = SoilTemperature_to10cm, y = MaxPercentIdentity)) +
  geom_point() + 
  facet_wrap(~Gene)

spearman_maxpident_logtemp <- stats.tidy %>% group_by(Gene) %>% do(tidy(cor.test(.$AveragePercentIdentity, log(.$SoilTemperature_to10cm), method = "spearman")))

qual <- c("adeB", "ClassA", "ClassB", "ClassC", "dfra12", "intI", "rplB", "sul2", "tolC", "vanA", "vanH", "vanX", "vanZ")

lala <- subset(stats.tidy, Gene %in% qual) %>% group_by(Gene) %>% do(tidy(cor.test(.$AverageLength, log(.$SoilTemperature_to10cm), method = "spearman")))


