#RNA Scope analysis with results from qupath script

library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(tidyverse)



setwd("/fh/fast/cheung_k/AndieX/Ami_image_analysis/Angptl7_RNAscope_3.2_test")


#read in measurement data as a list of dataframes
measurement_files_rna <- list.files(pattern = "*measurements")
measurement_data_rna <- lapply(measurement_files_rna, function(x)read.table(x, header=T, sep = "\t"))
#combine into one dataframe
measurement.df_rna = bind_rows(measurement_data_rna)

#read in annotation data as a list of dataframes
annotation_files_rna <- list.files(pattern = "*annotations.txt")
annotation_data_rna <- lapply(annotation_files_rna, function(x)read.table(x, header=T, sep = "\t"))
#combine into one dataframe
annotation.df_rna = bind_rows(annotation_data_rna)


#create rat id column
measurement.df_rna <- measurement.df_rna %>% 
  mutate(Rat = substr(Image,1,4))

annotation.df_rna <- annotation.df_rna %>% 
  mutate(Rat = substr(Image,1,4))

#add condition data
measurement.df_rna <- measurement.df_rna %>%
  mutate(Condition = case_when((Rat == 4427 | Rat == 4425 | Rat == 4426 | Rat == 4428) ~ 'Day 13', 
                               (Rat == 4423 | Rat == 4433 | Rat == 4431 | Rat == 4432) ~ 'Day 17',
                               (Rat == 4430 | Rat == 4437 | Rat == 4424 | Rat == 4422) ~ 'Day 22',
                               (Rat == 4435 | Rat == 4434 | Rat == 4429 | Rat == 4436) ~ 'Day 27'))

annotation.df_rna <- annotation.df_rna %>%
  mutate(Condition = case_when((Rat == 4427 | Rat == 4425 | Rat == 4426 | Rat == 4428) ~ 'Day 13', 
                               (Rat == 4423 | Rat == 4433 | Rat == 4431 | Rat == 4432) ~ 'Day 17',
                               (Rat == 4430 | Rat == 4437 | Rat == 4424 | Rat == 4422) ~ 'Day 22',
                               (Rat == 4435 | Rat == 4434 | Rat == 4429 | Rat == 4436) ~ 'Day 27'))


########

#filter for positive cells
pos_measurement.df <- measurement.df_rna %>%
  filter(Class == "Positive")


#plot distance to tumor boundary - not normalized (see below)
ggplot(pos_measurement.df, aes(x = Condition, y = Distance.to.annotation.Tumor.Inverse.µm))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.1, binwidth=20)


#set comparisons for stats
my_comparisons <- list( c("Day 13", "Day 17"), c("Day 17", "Day 22"), c("Day 22", "Day 27"), c("Day 13", "Day 27"))

ggplot(pos_measurement.df, aes(x = Condition, y = Distance.to.annotation.Tumor.Inverse.µm, color = Condition))+
  geom_violin(size = 1, scale = "width")+
  theme_bw()+ 
  geom_boxplot(width = 0.2, position =  position_dodge(width = 0.9), outlier.shape = NA, fill = "NA", size = .8) +
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Distance to Tumor Boundary (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons)


#plot distance to necrosis - not normalized (see below)
ggplot(pos_measurement.df, aes(x = Condition, y = Distance.to.annotation.Necrosis.µm))+
  geom_violin()

ggplot(pos_measurement.df, aes(x = Condition, y = Distance.to.annotation.Necrosis.µm, color = Condition))+
  geom_violin(size = 1, alpha = 0.5, scale = "width")+
  theme_bw()+ 
  geom_boxplot(width = 0.15, position =  position_dodge(width = 0.9), outlier.shape = NA, fill = "NA", size = .8) +
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Distance to Nearest Necrosis Boundary (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons)




#######

#overlay distance to necrosis and distance to tumor boundary

#first pivot data
pos_measurement.df.simple <- pos_measurement.df %>%
  select(c(Condition, Rat, Distance.to.annotation.Necrosis.µm, Distance.to.annotation.Tumor.Inverse.µm))

pivot_pos <- pos_measurement.df.simple %>%
  pivot_longer(cols = c(Distance.to.annotation.Necrosis.µm, Distance.to.annotation.Tumor.Inverse.µm),
               names_to = "Distance_to", values_to = "Microns")

#plot
ggplot(pivot_pos, aes(x= Condition, y=  Microns, color = Distance_to))+
  theme_bw()+
  geom_violin(scale = "width", position = position_dodge(width = 0.9), size = 1) +
  scale_color_manual(labels = c("Distance to Necrosis", "Distance to Tumor Edge"), values = c("#13783D", "#882256"), name = "Distance To:") +
  geom_boxplot(width = 0.2, position =  position_dodge(width = 0.9), outlier.shape = NA, size = 0.8) +
  ylab("Distance (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank())



########

#get summary data
pos_measurement.df.summary <- pos_measurement.df %>%
  group_by(Condition, Rat) %>%
  summarise(n_pos = n(), avg_tumor_dist = mean(Distance.to.annotation.Tumor.Inverse.µm),
            avg_necrosis_dist = mean(Distance.to.annotation.Necrosis.µm))

annotation.df.summary_rna <- annotation.df_rna %>%
  group_by(Condition, Rat) %>%
  summarise(total_necrotic_area = sum(Area.µm.2[Class=="Necrosis"]), tumor_area = sum(Area.µm.2[Class=="Tumor"]))

#combine annotation and measurement summaries
summary_data_rna <- left_join(pos_measurement.df.summary, annotation.df.summary_rna)

#plot with mean points
ggplot(summary_data_rna, aes(Condition, total_necrotic_area, color = Condition))+
  geom_violin(size = 1, alpha = 0.5, scale = "width")+
  theme_bw()+ 
  geom_point()+
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Total Necrotic Area (µm^2)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")


ggplot(summary_data_rna, aes(Condition, n_pos, color = Condition))+
  geom_violin(size = 1, alpha = 0.5, scale = "width")+
  theme_bw()+ 
  geom_point()+
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Number of Angptl7+ Cells")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test")

#plot distances summaries individually (not normalized)
ggplot(summary_data_rna, aes(Condition, avg_necrosis_dist, color = Condition))+
  geom_violin(size = 1, alpha = 0.5, scale = "width")+
  theme_bw()+ 
  geom_point()+
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Average Distance to Nearest Necrosis Boundary (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 17), axis.title.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons)

ggplot(summary_data_rna, aes(Condition, avg_tumor_dist, color = Condition))+
  geom_violin(size = 1, alpha = 0.5, scale = "width")+
  theme_bw()+ 
  geom_point()+
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Average Distance to Tumor Boundary (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 17), axis.title.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons)



#compare summary of necrosis to number of Angptl7+ cells
ggplot(summary_data_rna, aes(total_necrotic_area, n_pos))+
  geom_point() +
  theme_bw()+
  geom_smooth(method = "lm", alpha = .15, color = "#882256") +
  ylab("Number of Angptl+ Cells") + xlab("Total Necrotic Area (µm^2)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))


#get r square
eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

eq(summary_data_rna$total_necrotic_area, summary_data_rna$n_pos)

#####combined violin with summary (not normalized)

#first pivot data
summary_data.simple <- summary_data_rna %>%
  select(c(Condition, Rat, avg_tumor_dist, avg_necrosis_dist))

pivot_summary <- summary_data.simple %>%
  pivot_longer(cols = c(avg_tumor_dist, avg_necrosis_dist),
               names_to = "Distance_to", values_to = "Microns")

#plot
ggplot(pivot_summary, aes(x= Condition, y=  Microns, color = Distance_to))+
  theme_bw()+
  geom_violin(scale = "width", position = position_dodge(width = 1)) +
  scale_color_manual(labels = c("Distance to Necrosis", "Distance to Tumor Edge"), values = c("#13783D", "#882256"), name = "Distance To:") +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1)) +
  ylab("Distance (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank()) 





#######

#figure out percent of total cells that are positive
#first count positive and negative cells
measurement.df_perc <- measurement.df_rna %>%
  group_by(Condition, Rat, Class) %>%
  summarise(n_cells = n())

#divide positive cells by sum of positive and negative cells
measurement.df_perc <- measurement.df_perc %>%
  group_by(Condition, Rat) %>%
  summarise(percent_pos = n_cells[Class == "Positive"] / ((n_cells[Class == "Negative"]) + (n_cells[Class == "Positive"])))

#add to summary data
summary_data_rna <- left_join(summary_data_rna, measurement.df_perc)


#plot
ggplot(summary_data_rna, aes(Condition, percent_pos, color = Condition))+
  geom_violin(size = 1, alpha = 0.5, scale = "width")+
  theme_bw()+ 
  geom_point() +
  scale_color_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) +
  ylab("Percent Angptl7 Cells")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank()) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") 




#compare summary of necrosis to percent Angptl7+ cells
ggplot(summary_data_rna, aes(total_necrotic_area, percent_pos))+
  geom_point() +
  theme_bw()+
  geom_smooth(method = "lm", alpha = .15, color = "#882256") +
  ylab("Percent Angptl+ Cells") + xlab("Total Necrotic Area (µm^2)") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))

#r square value
eq(summary_data_rna$total_necrotic_area, summary_data_rna$percent_pos)

write.csv(summary_data, file = "Angptl7_summary_data_andie_analysis.csv")


########

#sample 15,000 random cells from each image and measure distance

measurement_random <- measurement.df_rna %>%
  group_by(Rat) %>%
  slice_sample(n=15000)


#overlay distance to necrosis and distance to tumor boundary
#first pivot data
measurement_random_simple <- measurement_random %>%
  select(c(Condition, Rat, Distance.to.annotation.Necrosis.µm, Distance.to.annotation.Tumor.Inverse.µm))

pivot_random <- measurement_random_simple %>%
  pivot_longer(cols = c(Distance.to.annotation.Necrosis.µm, Distance.to.annotation.Tumor.Inverse.µm),
               names_to = "Distance_to", values_to = "Microns")

#plot
ggplot(pivot_random, aes(x= Condition, y=  Microns, color = Distance_to))+
  theme_bw()+
  geom_violin(scale = "width", position = position_dodge(width = 0.9), size = 1) +
  scale_color_manual(labels = c("Distance to Necrosis", "Distance to Tumor Edge"), values = c("#13783D", "#882256"), name = "Distance To:") +
  geom_boxplot(width = 0.2, position =  position_dodge(width = 0.9), outlier.shape = NA, size = 0.8) +
  ylab("Distance (µm)")+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank())


ggplot(measurement_random_simple, aes(x = Condition, y = Distance.to.annotation.Necrosis.µm))+
  geom_violin()



###################

### Normalizing to null density distribution

#make list of images
Images <- unique(measurement_random$Rat)

#remove rat with no necrosis
Images <- Images[Images != "4428"]

#loop through each and calculate density and then make a function for it
out_functions <- vector('list')

for (i in Images){

  df <- measurement_random %>%
    filter((Rat == i) & (!is.na(Distance.to.annotation.Necrosis.µm)))
  
  min_value <- min(df$Distance.to.annotation.Necrosis.µm, na.rm = TRUE)
  max_value <- max(df$Distance.to.annotation.Necrosis.µm, na.rm = TRUE)
  count_value <- nrow(df)

  kde <- density(x = df$Distance.to.annotation.Necrosis.µm, from= min_value, to = max_value,  n = count_value, na.rm = TRUE)

  plot(kde)
  
  kde_fun <- approxfun(kde$x, kde$y)
  
  out_functions[[i]] <- kde_fun
}


#loop to add null density at each distance in positive data detections, must split by image
pos_measurement_kde_list <- vector('list')

for (i in Images){
  
  kde_fun <- out_functions[[i]]
  
  df <- pos_measurement.df %>%
    filter((Rat == i) & (!is.na(Distance.to.annotation.Necrosis.µm)))
  
  df <- df %>%
    mutate(kde_dist = kde_fun(Distance.to.annotation.Necrosis.µm))
  
  pos_measurement_kde_list[[i]] <- df
}

#make back into one dataframe
kde_pos = bind_rows(pos_measurement_kde_list)

#remove na
kde_pos <- kde_pos %>%
  filter(!is.na(kde_dist))

#summarize results
kde_pos_summary <- kde_pos %>%
  group_by(Condition, Rat) %>%
  summarize(weighted_avg_distance_necrosis = weighted.mean(Distance.to.annotation.Necrosis.µm, kde_dist),
         avg_distance_necrosis = mean(Distance.to.annotation.Necrosis.µm))

#check
ggplot(kde_pos_summary, aes(x = Condition, y = weighted_avg_distance_necrosis))+
  geom_violin() + geom_point()

ggplot(kde_pos_summary, aes(x = Condition, y = avg_distance_necrosis))+
  geom_violin() + geom_point()


#compare to random 
random_measurement_kde_list <- vector('list')

for (i in Images){
  
  kde_fun <- out_functions[[i]]
  
  df <- measurement_random %>%
    filter((Rat == i) & (!is.na(Distance.to.annotation.Necrosis.µm)))
  
  df <- df %>%
    mutate(kde_dist = kde_fun(Distance.to.annotation.Necrosis.µm))
  
  random_measurement_kde_list[[i]] <- df
}

#make back into one dataframe
kde_random = bind_rows(random_measurement_kde_list)

#remove na
kde_random <- kde_random %>%
  filter(!is.na(kde_dist))

#summarize results
kde_random_summary <- kde_random %>%
  group_by(Condition, Rat) %>%
  summarize(weighted_avg_distance_necrosis = weighted.mean(Distance.to.annotation.Necrosis.µm, kde_dist),
            avg_distance_necrosis = mean(Distance.to.annotation.Necrosis.µm))

#check
ggplot(kde_random_summary, aes(x = Condition, y = weighted_avg_distance_necrosis))+
  geom_violin() + geom_point()

ggplot(kde_random_summary, aes(x = Condition, y = avg_distance_necrosis))+
  geom_violin() + geom_point()


ggplot(kde_random, aes (x = Rat, y = kde_dist)) +
  geom_violin()


#### do the same with distance to tumor


#make list of images
Images <- unique(measurement_random$Rat)

#remove rat with no necrosis
Images <- Images[Images != "4428"]

#loop through each and calculate density and then make a function for it
out_functions_tumor <- vector('list')

for (i in Images){
  
  df <- measurement_random %>%
    filter((Rat == i) & (!is.na(Distance.to.annotation.Tumor.Inverse.µm)))
  
  min_value <- min(df$Distance.to.annotation.Tumor.Inverse.µm, na.rm = TRUE)
  max_value <- max(df$Distance.to.annotation.Tumor.Inverse.µm, na.rm = TRUE)
  count_value <- nrow(df)
  
  kde <- density(x = df$Distance.to.annotation.Tumor.Inverse.µm, from= min_value, to = max_value,  n = count_value, na.rm = TRUE)
  
  plot(kde)
  
  kde_fun <- approxfun(kde$x, kde$y)
  
  out_functions_tumor[[i]] <- kde_fun
}


#loop to add null density at each distance in positive data detections, must split by image
pos_measurement_kde_list_tumor <- vector('list')

for (i in Images){
  
  kde_fun <- out_functions_tumor[[i]]
  
  df <- pos_measurement.df %>%
    filter((Rat == i) & (!is.na(Distance.to.annotation.Tumor.Inverse.µm)))
  
  df <- df %>%
    mutate(kde_dist = kde_fun(Distance.to.annotation.Tumor.Inverse.µm))
  
  pos_measurement_kde_list_tumor[[i]] <- df
}

#make back into one dataframe
kde_pos_tumor = bind_rows(pos_measurement_kde_list_tumor)

#remove na kde_dist
kde_pos_tumor <- kde_pos_tumor %>%
  filter(!is.na(kde_dist))

#summarize results
kde_pos_summary_tumor <- kde_pos_tumor %>%
  group_by(Condition, Rat) %>%
  summarize(weighted_avg_distance_tumor = weighted.mean(Distance.to.annotation.Tumor.Inverse.µm, kde_dist),
            avg_distance_tumor = mean(Distance.to.annotation.Tumor.Inverse.µm))

#check
ggplot(kde_pos_summary_tumor, aes(x = Condition, y = weighted_avg_distance_tumor))+
  geom_violin() + geom_point()

ggplot(kde_pos_summary_tumor, aes(x = Condition, y = avg_distance_tumor))+
  geom_violin() + geom_point()


#compare to random 
random_measurement_kde_list_tumor <- vector('list')

for (i in Images){
  
  kde_fun <- out_functions_tumor[[i]]
  
  df <- measurement_random %>%
    filter((Rat == i) & (!is.na(Distance.to.annotation.Tumor.Inverse.µm)))
  
  df <- df %>%
    mutate(kde_dist = kde_fun(Distance.to.annotation.Tumor.Inverse.µm))
  
  random_measurement_kde_list_tumor[[i]] <- df
}

#make back into one dataframe
kde_random_tumor = bind_rows(random_measurement_kde_list)

#remove na kde_dist
kde_random_tumor <- kde_random_tumor %>%
  filter(!is.na(kde_dist))

#summarize results
kde_random_summary_tumor <- kde_random_tumor %>%
  group_by(Condition, Rat) %>%
  summarize(weighted_avg_distance_tumor = weighted.mean(Distance.to.annotation.Tumor.Inverse.µm, kde_dist),
            avg_distance_tumor = mean(Distance.to.annotation.Tumor.Inverse.µm))

#check
ggplot(kde_random_summary_tumor, aes(x = Condition, y = weighted_avg_distance_tumor))+
  geom_violin() + geom_point()

ggplot(kde_random_summary_tumor, aes(x = Condition, y = avg_distance_tumor))+
  geom_violin() + geom_point()







#######  account for changes in size by subtracting random from positive

kde_sum_tumor <- left_join(kde_random_summary_tumor, kde_pos_summary_tumor, by = c("Condition", "Rat"), suffix = c(".null", ".pos"))

kde_sum_tumor <- kde_sum_tumor %>%
  mutate(pos_minus_null.tumor = weighted_avg_distance_tumor.pos- weighted_avg_distance_tumor.null)

ggplot(kde_sum_tumor, aes(x = Condition, y = pos_minus_null.tumor))+
  geom_violin() + geom_point()


#do the same with necrosis

kde_sum_necrosis <- left_join(kde_random_summary, kde_pos_summary, by = c("Condition", "Rat"), suffix = c(".null", ".pos"))

kde_sum_necrosis <- kde_sum_necrosis %>%
  mutate(pos_minus_null.necrosis = weighted_avg_distance_necrosis.pos- weighted_avg_distance_necrosis.null)

ggplot(kde_sum_necrosis, aes(x = Condition, y = pos_minus_null.necrosis))+
  geom_violin() + geom_point()




#plot distance to necrosis and tumor together 
#combine dataframes
kde_sum_both <- left_join(kde_sum_tumor, kde_sum_necrosis, by = c("Condition", "Rat"))

#isolate data you want
kde_sum_both_simple <- kde_sum_both %>%
  select(c(Condition, Rat, pos_minus_null.tumor, pos_minus_null.necrosis))


#pivot
pivot_kde_summary <- kde_sum_both_simple %>%
  pivot_longer(cols = c(pos_minus_null.tumor, pos_minus_null.necrosis),
               names_to = "Distance_to", values_to = "Normalized_dist")


#figure out getting se 

#plot
ggplot(pivot_kde_summary, aes(x = Condition, y = Normalized_dist, fill = Distance_to)) +
  geom_violin()

#average
sum_kde_pivot <- pivot_kde_summary %>%
  group_by(Condition, Distance_to)%>%
  summarize(mean_dist = mean(Normalized_dist), sd = sd(Normalized_dist), n= n(), se = sd/sqrt(n))



#plot
ggplot(sum_kde_pivot, aes(x = Condition, y = mean_dist, fill = Distance_to)) +
  geom_col(position = "Dodge") + theme_bw() +
  ylab("Average distance away compared to null (microns)") +
  scale_fill_manual(labels = c("Distance to Necrosis", "Distance to Tumor Edge"), values = c("#13783D", "#882256"), name = "Distance To:") +
  geom_errorbar(aes(ymin = mean_dist-se, ymax = mean_dist+se), position = "Dodge")+
  geom_text(aes(y = 550, label = round(mean_dist, digits = 1)), position = position_dodge(0.9), size = 2.5)
  
  

############################

#make enrichment score so that closest to necrosis = 1, closest to tumor border  = -1

kde_sum_both_es <- kde_sum_both %>%
  mutate(enrichment_score = ifelse(pos_minus_null.necrosis <= 0, ((weighted_avg_distance_necrosis.null - weighted_avg_distance_necrosis.pos)/weighted_avg_distance_necrosis.null), 
                                   ((weighted_avg_distance_tumor.pos - weighted_avg_distance_tumor.null)/weighted_avg_distance_tumor.null)))


#summarize with stats
kde_sum_both_es_stats <- kde_sum_both_es %>%
  group_by(Condition)%>%
  summarize(mean_es = mean(enrichment_score), sd = sd(enrichment_score), n= n(), se = sd/sqrt(n))

kde_sum_both_es<- kde_sum_both_es %>%
  mutate(pos_minus_mutliplied = pos_minus_null.necrosis* pos_minus_null.tumor)

#plot
ggplot(kde_sum_both_es_stats, aes(x = Condition, y = mean_es, fill = Condition)) +
  geom_col() + theme_bw() +
  ylab("Necrotic Enrichment Score") +
  scale_fill_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) + 
  geom_errorbar(aes(ymin = mean_es-se, ymax = mean_es+se))+
  geom_text(aes(y = 0.75, label = round(mean_es, digits = 1)), position = position_dodge(0.9), size = 2.5)+
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank())


ggplot(kde_sum_both_es, aes(x = Condition, y = enrichment_score, fill = Condition)) +
  geom_violin() + geom_point() + theme_bw() + 
  ylab("Necrotic Enrichment Score") +
  scale_fill_manual(values = c("#332F85", "#13783D", "#44AA99", "#882256")) + 
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18), axis.title.x = element_blank())+
  stat_compare_means(comparisons = my_comparisons) 






