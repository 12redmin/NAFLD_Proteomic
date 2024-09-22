

library(tidyverse)
library(feather)
library(dplyr)
library(ggplot2)
library(stringr)
library(arrow)
library(tableone)
library(RANN)


newdata <- read_feather("F:/UKB_newdata/ukb675591.feather")

olink_data<- read.table("E:/UKBPPP/Omics_data/olink_data.txt",sep = "\t",header = T)

message("N eids: ",olink_data %>% distinct(eid) %>% nrow)
message("N proteins: ",olink_data %>% distinct(protein_id) %>% nrow)
#N eids: 53062; N proteins: 2923

#only select instance=0
olink_data<-olink_data[olink_data$ins_index==0,]
message("N eids: ",olink_data %>% distinct(eid) %>% nrow)
message("N proteins: ",olink_data %>% distinct(protein_id) %>% nrow)
#N eids: 53018; N proteins: 2923
olink_eid<-unique(olink_data$eid)
olink_eid<-data.frame(unique(olink_data$eid))
names(olink_eid)[1] <- "eid" 



protein_ids = read.csv("E:/UKBPPP/Omics_data/UKBPPP/protein_id_code.csv", header=T) %>%
  tidyr::separate(protein, sep = ";", into = c("protein","full_name")) 
write.csv(protein_ids,"E:/UKBPPP/Omics_data/UKBPPP/protein_detail_name.csv")


olink_data = olink_data %>%
  filter(protein_id %in% protein_ids$protein_id) %>%
  left_join(protein_ids,by="protein_id")

dat_long <- readRDS("E:/UKBPPP/Omics_data/UKBPPP/proteomics__instance0.rds")



dat_wide <- olink_data %>%
  dplyr::select(eid,protein,result) %>%
  pivot_wider(id_cols = eid, values_from = result, names_from = protein)


olink_NA <- data.frame(
  Protein = character(0),
  Missing_Percentage = numeric(0)
)


for (col_name in colnames(dat_wide)[-1]) {  # 排除第一列（ID列）
  
  missing_percentage <- mean(is.na(dat_wide[[col_name]]))
  
  olink_NA <- rbind(olink_NA, data.frame(protein = col_name, Missing_Percentage = missing_percentage))
}


olink_NA$keep<-ifelse(olink_NA$Missing_Percentage>0.2,0,1)

olink_instance0<-merge(protein_ids,olink_NA,by="protein")

table(olink_instance0$keep==0)
#FALSE  TRUE 
#2911    12 


data_baseline<-readRDS("E:/UKBPPP/Omics_data/UKBPPP/data_baseline.rds")

dat_wide <- readRDS("E:/UKBPPP/Omics_data/UKBPPP/proteomics_wide_instance0.rds")

olink_instance0<-read.csv("E:/UKBPPP/Omics_data/UKBPPP/olink_instance0.csv",header=T)


olink_remove<-olink_instance0[olink_instance0$keep==0,]

remove_cols <- intersect(names(dat_wide), olink_remove$protein)

dat_wide <- dat_wide[, !names(dat_wide) %in% remove_cols]
rm(olink_remove,remove_cols)

dat_wide<-subset(dat_wide,dat_wide$eid%in%data_baseline$eid,)


data_imputation_3 <- readRDS("E:/UKBPPP/Omics_data/UKBPPP/data_imputation_3.rds")



data_demographic<-data_imputation_3 %>% dplyr::select(eid,
                                                      age,
                                                      sex,
                                                      region,
                                                      tdi,
                                                      ethnic,
                                                      education)


data_demographic_female<-data_demographic[data_demographic$sex=="0",]

data_demographic_male<-data_demographic[data_demographic$sex=="1",]

dat_wide_female<-dat_wide[dat_wide$eid%in%data_demographic_female$eid,]

dat_wide_male<-dat_wide[dat_wide$eid%in%data_demographic_male$eid,]




calculate_protein_correlation <- function(data_protein) {
  protein_data <- data_protein[, -1]
  protein_data <- na.omit(protein_data)
  protein_correlation <- cor(protein_data)
  return(protein_correlation)
}



protein_correlation_female <- calculate_protein_correlation(dat_wide_female)

protein_correlation_male <- calculate_protein_correlation(dat_wide_male)

saveRDS(protein_correlation_female,"E:/UKBPPP/Omics_data/UKBPPP/protein_correlation_female.rds",compress = F)
saveRDS(protein_correlation_male,"E:/UKBPPP/Omics_data/UKBPPP/protein_correlation_male.rds",compress = F)
#rm(protein_correlation_female,protein_correlation_male)


select_related_proteins <- function(protein_correlation, n = 10) {
  selected_proteins <- list()
  for (i in 1:ncol(protein_correlation)) {
    correlated_proteins <- sort(protein_correlation[, i], decreasing = TRUE, na.last = TRUE)
    correlated_proteins <- correlated_proteins[!is.na(correlated_proteins)]
    selected_proteins[[i]] <- names(correlated_proteins)[-i][1:min(n, length(correlated_proteins)-1)]
  }
  return(selected_proteins)
}


selected_proteins_female <- select_related_proteins(protein_correlation_female)

selected_proteins_male <- select_related_proteins(protein_correlation_male)




imputed_data_female <- dat_wide_female[,1]

for(i in 1:2912){
  protein_imp <- dat_wide_female[, c(1, i)]
  selected_proteins_imp <- dat_wide_female[c("eid", selected_proteins_female[[i-1]])]
  combine_data_imp <- merge(data_demographic_female, protein_imp, by = "eid")
  combine_data_imp <- merge(combine_data_imp, selected_proteins_imp, by = "eid")
  
  protein_imp_female <- knnImputation(combine_data_imp, k = 10, scale = TRUE, meth = "weighAvg", distData = NULL)
  
  
  protein_imp_female <- protein_imp_female[,c(1,8)]
  
  imputed_data_female <-merge(imputed_data_female,protein_imp_female,by="eid")
  
}


for(i in 1:2912){
  protein_imp <- dat_wide_male[, c(1, i)]
  selected_proteins_imp <- dat_wide_male[c("eid", selected_proteins_male[[i-1]])]
  combine_data_imp <- merge(data_demographic_male, protein_imp, by = "eid")
  combine_data_imp <- merge(combine_data_imp, selected_proteins_imp, by = "eid")
  protein_imp_male <- knnImputation(combine_data_imp, k = 10, scale = TRUE, meth = "weighAvg", distData = NULL)
  protein_imp_male <- protein_imp_male[,c(1,8)]
  imputed_data_male <-merge(imputed_data_male,protein_imp_male,by="eid")
  
}

dat_wide_imputed<-merge(dat_wide[,1],rbind(imputed_data_female,imputed_data_male),by="eid")

saveRDS(dat_wide_imputed,"E:/UKBPPP/Omics_data/UKBPPP/dat_wide_imputed.rds",compress = F)