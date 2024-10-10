

dat_wide_imputed<-readRDS("E:/UKBPPP/Omics_data/UKBPPP/dat_wide_imputed.rds")
data_imputation_3 <- readRDS("E:/UKBPPP/Omics_data/UKBPPP/data_imputation_3.rds")


#Cox----------------------------------
NAFLD_cox_model1 <- data.frame(
  Protein_Name = character(),
  HR = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)
NAFLD_cox_model2 <- data.frame(
  Protein_Name = character(),
  HR = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

NAFLD_cox_model3 <- data.frame(
  Protein_Name = character(),
  HR = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)





for (i in 30:2940) { 
  protein_name <- names(data_model)[i]
  
  fit_1 <- coxph(Surv(NAFLDtime,NAFLD) ~ age + sex+ data_model[, i], data = data_model,id=eid)
  
  
  summary_fit <- summary(fit_1)
  coef <- summary_fit$coefficients["data_model[, i]", "coef"]
  HR <- round(exp(coef), digits = 2)
  CI <- round(exp(confint(fit_1, "data_model[, i]")), digits = 2)
  P_Value <- summary_fit$coefficients["data_model[, i]", "Pr(>|z|)"]
  

  NAFLD_cox_model1 <- rbind(NAFLD_cox_model1, data.frame(
    Protein_Name = protein_name,
    HR = HR,
    CI_Lower = CI[1],
    CI_Upper = CI[2],
    P_Value = P_Value
  ))
  
  
  fit_2 <- coxph(Surv(NAFLDtime,NAFLD) ~ age + sex + ethnic +  tdi + education + smoke + drink + bmi + PA + 
                   hypertension + diabetes + data_model[, i], data = data_model,id=eid)
  
  
  NAFLD_cox_model2 <- rbind(NAFLD_cox_model2, data.frame(
    Protein_Name = protein_name,
    HR = HR,
    CI_Lower = CI[1],
    CI_Upper = CI[2],
    P_Value = P_Value
  ))
  
  
  fit_3 <- coxph(Surv(NAFLDtime,NAFLD) ~ age + sex + ethnic +  tdi + education + smoke + drink + bmi + PA + 
                   hypertension + diabetes + ALP + ALT + AST + GGT+ data_model[, i], data = data_model,id=eid)
  
  
  NAFLD_cox_model3 <- rbind(NAFLD_cox_model3, data.frame(
    Protein_Name = protein_name,
    HR = HR,
    CI_Lower = CI[1],
    CI_Upper = CI[2],
    P_Value = P_Value
  ))
  
  
  
}

NAFLD_cox_model1$HR_CI<-paste0(NAFLD_cox_model1$HR," (",NAFLD_cox_model1$CI_Lower,"-",NAFLD_cox_model1$CI_Upper,")")
NAFLD_cox_model1$Pvalue_adjust<-p.adjust(NAFLD_cox_model1$P_Value,method = "bonferroni")
table(NAFLD_cox_model1$Pvalue_adjust<0.05)


NAFLD_cox_model2$HR_CI<-paste0(NAFLD_cox_model2$HR," (",NAFLD_cox_model2$CI_Lower,"-",NAFLD_cox_model2$CI_Upper,")")
NAFLD_cox_model2$Pvalue_adjust<-p.adjust(NAFLD_cox_model2$P_Value,method = "bonferroni")
table(NAFLD_cox_model2$Pvalue_adjust<0.05)

NAFLD_cox_model3$HR_CI<-paste0(NAFLD_cox_model3$HR," (",NAFLD_cox_model3$CI_Lower,"-",NAFLD_cox_model3$CI_Upper,")")
NAFLD_cox_model3$Pvalue_adjust<-p.adjust(NAFLD_cox_model3$P_Value,method = "bonferroni")
table(NAFLD_cox_model3$Pvalue_adjust<0.05)


NAFLD_cox_sig_1<-NAFLD_cox_model1[NAFLD_cox_model1$Pvalue_adjust<0.05,]
NAFLD_cox_sig_2<-NAFLD_cox_model2[NAFLD_cox_model2$Pvalue_adjust<0.05,]
NAFLD_cox_sig_3<-NAFLD_cox_model3[NAFLD_cox_model3$Pvalue_adjust<0.05,]


dat_wide_shared_2_final<-subset(dat_wide_imputed,,names(dat_wide_imputed)%in%NAFLD_cox_sig_1$Protein_Name
                                &names(dat_wide_imputed)%in%NAFLD_cox_sig_2$Protein_Name
                                &names(dat_wide_imputed)%in%NAFLD_cox_sig_3$Protein_Name
                                |names(dat_wide_imputed)=="eid")

NAFLD_outcome<-subset(data_imputation_3,,c(eid,NAFLD))
dat_sig_shared_2_final<-merge(NAFLD_outcome,dat_wide_shared_2_final,by="eid")
dat_sig_shared_2_final$eid<-NULL

saveRDS(dat_sig_shared_2_final,"E:/UKBPPP/Final_result/NAFLD_dat_sig_shared_2_final.rds")


library(reportROC)
result_ROC_2 <- data.frame()


for (i in 2:325) { 
  protein_name <- names(dat_sig_shared_2_final)[i]
  
  ROC.info <- reportROC(gold=dat_sig_shared_2_final$NAFLD,
                        predictor=dat_sig_shared_2_final[,i],
                        important="se",
                        plot=TRUE)
  
  
  ROC.info$Protein_Name<-protein_name
  
  
  result_ROC_2 <- rbind(result_ROC_2,ROC.info)
  
}


NAFLD_result_ROC_2<-subset(result_ROC_2,,c(Protein_Name,1:39))
NAFLD_result_ROC_2$outcome<-"NAFLD"
NAFLD_result_ROC_2<-subset(NAFLD_result_ROC_2,,c(Protein_Name,outcome,2:40))

result_ROC_2$outcome<-"NAFLD"
write.csv(result_ROC_2,"E:/UKBPPP/Final_result/NAFLD_result_ROC_2.csv")

table(result_ROC_2$AUC>=0.6)



#EN-----------------------------------------------------
library(glmnet)
library(foreign)
library(RColorBrewer)

dat_sig_shared_2_final<-readRDS("E:/UKBPPP/Final_result/NAFLD_dat_sig_shared_2_final.rds")


NAFLD_result_ROC_2<-read.csv("E:/UKBPPP/Final_result/NAFLD_result_ROC_2.csv")

NAFLD_result_ROC_2<-subset(NAFLD_result_ROC_2,NAFLD_result_ROC_2$AUC>0.6)
#271
table(NAFLD_result_ROC_2$AUC>=0.6)

dat_sig_shared_2_final<-subset(dat_sig_shared_2_final,,names(dat_sig_shared_2_final)=="NAFLD"|
                                 names(dat_sig_shared_2_final)%in%NAFLD_result_ROC_2$Protein_Name)



x2 <- as.matrix(dat_sig_shared_2_final[, 2:272])
x2 <- scale(x2)

dat_sig_shared_2_final$NAFLD <- as.numeric(dat_sig_shared_2_final$NAFLD)
y2<-as.matrix(dat_sig_shared_2_final[,1])



model_lasso2 <- glmnet(x2, y2, family="binomial", nlambda=100, alpha=0.05)
print(model_lasso2)


colors <- brewer.pal(5, "Set1")
plot(model_lasso2, xvar = "lambda", label = FALSE, lwd = 1.5, col = colors)



set.seed(123)
cvfit2 <- cv.glmnet(x2, y2,family = "binomial",alpha = 0.5)
plot(cvfit2)



lambda.min2 <- cvfit2$lambda.min
lambda.se2 <- cvfit2$lambda.1se
lambda.min2
#0.001094151
lambda.se2
#0.01119897



model_lasso_min2 <- glmnet(x2, y2,family = "binomial",alpha = 0.5, lambda = lambda.min2,exact = F)
model_lasso_se2 <- glmnet(x2, y2,family = "binomial",alpha = 0.5, lambda = lambda.se2,exact = F)



model_lasso_min2$beta
model_lasso_se2$beta

protein_min2 <- rownames(model_lasso_min2$beta)[as.numeric(model_lasso_min2$beta)!=0]
protein_se2 <- rownames(model_lasso_se2$beta)[as.numeric(model_lasso_se2$beta)!=0]
length(protein_min2)#67 proteins


length(protein_se2)#18 proteins

protein_min2
# [1] "TNFRSF10A" "TNFRSF11A" "TNFRSF14"  "HAO1"      "GZMA"      "GUSB"      "GSTA1"     "FABP4"     "F7"        "ERBB2"     "GGH"      
#[12] "GGT1"      "FGF21"     "NPY"       "NOMO1"     "INHBC"     "ITGA5"     "IL1RN"     "IL6"       "LEP"       "BST2"      "ARSA"     
#[23] "ASGR1"     "ACY1"      "ADAM8"     "ACE2"      "ADGRG1"    "ADH4"      "AGRN"      "CTSO"      "CXCL10"    "CREG1"     "CRLF1"    
#[34] "DCXR"      "CDHR2"     "CES1"      "CDCP1"     "CPM"       "CHI3L1"    "RET"       "RBKS"      "SIGLEC7"   "SEMA3F"    "PLIN1"    
#[45] "TGFBI"     "PRCP"      "INSL5"     "KRT8"      "LAMP2"     "TREH"      "RRM2"      "PRAP1"     "PSAP"      "SIGLEC8"   "SHISA5"   
#[56] "CEMIP2"    "CFI"       "ADAMTSL2"  "ADGRE1"    "FUOM"      "HJV"       "MAN2B2"    "ECHDC3"    "BPIFB2"    "FCAMR"     "TCOF1"    
#[67] "UPB1"  

protein_se2
# [1] "GUSB"     "IGSF3"    "GGT1"     "IL1RN"    "KRT18"    "ACY1"     "ACE2"     "ADGRG1"   "CES1"     "CPM"      "PALM2"    "PCBD1"   
#[13] "KRT8"     "PRAP1"    "ADAMTSL2" "FUOM"     "FTCD"     "IGSF9" 


NAFLD_ML_protein<-as.data.frame(protein_se2)
names(NAFLD_ML_protein)[1]<-"Protein"
NAFLD_ML_protein$outcome<-"NAFLD"




#LightGBM----------------------------------------------------------
library(bonsai)
library(lightgbm)
library(treesnip)
library(shapviz)
library(ggplot2)
library(tidymodels)
tidymodels_prefer()

NAFLD_protein_select<-readRDS("E:/UKBPPP/Final_result/NAFLD_protein_select.rds")

set.seed(123)
dat_split <- initial_split(NAFLD_protein_select, prop = 0.80)
dat_train <- training(dat_split)
dat_test  <-  testing(dat_split)





folds <- vfold_cv(dat_train,v=10)


lightgbm_rec <- recipe(NAFLD~.,data=dat_train)



lightgbm_tune_mod <-boost_tree(trees = tune(), tree_depth = tune(), 
                               learn_rate = tune(), min_n = tune(), loss_reduction = tune()) %>%
  set_engine('lightgbm') %>%
  set_mode('classification')


lightgbm_tune_wf <- workflow() %>% 
  add_model(lightgbm_tune_mod) %>% 
  add_recipe(lightgbm_rec)


lightgbm_grid <- grid_random(trees(),tree_depth(),learn_rate(), min_n(), loss_reduction(),size = 100)



lightgbm_tune_fit_NAFLD <- tune_grid(          
  lightgbm_tune_wf,          
  grid=lightgbm_grid,           
  resamples=folds)

#saveRDS(lightgbm_tune_fit_NAFLD,"E:/UKBPPP/Final_result/lightgbm_tune_fit_NAFLD.rds",compress = F)

lightgbm_tune_fit_NAFLD<-readRDS("E:/UKBPPP/Final_result/lightgbm_tune_fit_NAFLD.rds")

autoplot(lightgbm_tune_fit_NAFLD)

lightgbm_bst_NAFLD <- select_best(lightgbm_tune_fit_NAFLD,metric="roc_auc")

#REF---------------------------------------------------
for (i in 18:1) {
  
  if(i==1){
    selected_features <- importance_tune$Feature[1:i]
    
    
    ROC.info <- reportROC(gold=dat_train$NAFLD,
                          predictor=dat_train[, c(selected_features)],
                          important="se",
                          plot=TRUE)
    
    auc_value<-round(as.numeric(ROC.info$AUC), digits = 10)
    
    
    auc_lower<- ROC.info$AUC.low
    auc_upper<- ROC.info$AUC.up
    
    temp_ROC_df <- data.frame(selected_protein = importance_tune$Feature[i],
                              accumulate_auc =auc_value,
                              auc_Lower = auc_lower,
                              auc_Upper = auc_upper)
    
    
    result_accumulate_auc <- rbind(result_accumulate_auc, temp_ROC_df)
    
  }
  else{
    selected_features <- importance_tune$Feature[1:i]
    
    
    dat_protein_select <- NAFLD_protein_select[, c("NAFLD",selected_features)]
    
    set.seed(123)
    dat_split <- initial_split(dat_protein_select, prop = 0.80)
    dat_train <- training(dat_split)
    dat_test  <-  testing(dat_split)
    
    
    lgb_data <- as.matrix(dat_train)
    lgb_train <- lgb.Dataset(data = lgb_data[, -1], 
                             label = lgb_data[, 1])
    

    lightgbm_tune_fit_NAFLD<-readRDS("E:/UKBPPP/Final_result/lightgbm_tune_fit_NAFLD.rds")
    lightgbm_bst_NAFLD <- select_best(lightgbm_tune_fit_NAFLD,metric="roc_auc")
    lgb_params <- list(objective = 'binary',metric = 'auc',
                       n_estimators=lightgbm_bst_NAFLD$trees,
                       min_data_in_leaf=lightgbm_bst_NAFLD$min_n,
                       max_depth =lightgbm_bst_NAFLD$tree_depth,
                       learning_rate  = lightgbm_bst_NAFLD$learn_rate,
                       min_gain_to_split =lightgbm_bst_NAFLD$loss_reduction)   
    
    lgb_model<- lgb.train(
      params = lgb_params,
      data = lgb_train,
      valids = list(valid = lgb.Dataset(data = lgb_data[, -1], label = lgb_data[, 1]))
    )
    
    
    shap_lgb<- shapviz(lgb_model, X_pred = lgb_data[, -1])
    shap_values<-shap_lgb[["X"]]
    
    sv_importance(
      shap_lgb,
      kind = c("beeswarm"),
      max_display = 20L,
      fill = "#fca50a",
      bar_width = 2/3,
      bee_width = 0.4,
      bee_adjust = 0.5,
      viridis_args = getOption("shapviz.viridis_args"),
      color_bar_title = "Protein value",
      show_numbers = TRUE,
      format_fun = format_max,
      number_size = 3.2)
    
    
    lgb_predicted_prob_tr <- predict(lgb_model, as.matrix(dat_test[,-1]))
    roc <- roc(dat_test$NAFLD, lgb_predicted_prob_tr)
    auc_value <- as.numeric(roc$auc)
    ci_auc <- ci.auc(roc)
    
    
    temp_ROC_df <- data.frame(selected_protein = importance_tune$Feature[i],
                              accumulate_auc = ci_auc[2],
                              auc_Lower = ci_auc[1],
                              auc_Upper = ci_auc[3])
    
    result_accumulate_auc <- rbind(result_accumulate_auc, temp_ROC_df)
    
    
  } 
}
write.csv(result_accumulate_auc,"E:/UKBPPP/Final_result/NAFLD_result_accumulate_auc.csv")

result_accumulate_auc<-read.csv("E:/UKBPPP/Final_result/NAFLD_result_accumulate_auc.csv")


write.csv(shap_values,"E:/UKBPPP/Final_result/NAFLD_8protein_shap.csv")


importance_tune<-read.csv("E:/UKBPPP/Final_result/NAFLD_importance.csv")



NAFLD_protein_select<-readRDS("E:/UKBPPP/Final_result/NAFLD_protein_select.rds")

set.seed(123)
dat_split <- initial_split(NAFLD_protein_select, prop = 0.80)
dat_train <- training(dat_split)
dat_test  <-  testing(dat_split)

lgb_data <- as.matrix(dat_train)
lgb_train <- lgb.Dataset(data = lgb_data[, -1], 
                         label = lgb_data[, 1])
lightgbm_tune_fit_NAFLD<-readRDS("E:/UKBPPP/Final_result/lightgbm_tune_fit_NAFLD.rds")
lightgbm_bst_NAFLD <- select_best(lightgbm_tune_fit_NAFLD,metric="roc_auc")

lightgbm_bst_NAFLD_5<-show_best(lightgbm_tune_fit_NAFLD,metric="roc_auc")


delong_test <- data.frame(selected_protein = importance_tune$Feature[18],
                          Delong_P=0)

for (i in 18:2) {
  
  if(i==2){
    selected_features_1 <- importance_tune$Feature[1]
    
    roc1 <- roc(dat_train$NAFLD,dat_train[, c(selected_features_1)])
    
    
    
    selected_features_2 <- importance_tune$Feature[1:i]
    
    
    dat_protein_select <- as.matrix(dat_train[, c("NAFLD",selected_features_2)])
    
    dat_protein_select_train <- lgb.Dataset(data = dat_protein_select[, -1], 
                                            label = dat_protein_select[, 1])
    
    
    dat_protein_select_test<-dat_test[, c("NAFLD",selected_features_2)]
    
    
    # 参数    
    lgb_params <- list(objective = 'binary',metric = 'auc',
                       n_estimators=lightgbm_bst_NAFLD$trees,
                       min_data_in_leaf=lightgbm_bst_NAFLD$min_n,
                       max_depth =lightgbm_bst_NAFLD$tree_depth,
                       learning_rate  = lightgbm_bst_NAFLD$learn_rate,
                       min_gain_to_split =lightgbm_bst_NAFLD$loss_reduction)   
    
    
    lgb_model <- lgb.train(
      params = lgb_params,
      data = dat_protein_select_train,
      valids = list(valid = lgb.Dataset(data = dat_protein_select[, -1], label = dat_protein_select[, 1]))
    )
    
    
    lgb_predicted_prob_tr <- predict(lgb_model, as.matrix(dat_protein_select_test[,-1]))
    
    
    roc2 <- roc(dat_protein_select_test$NAFLD, lgb_predicted_prob_tr)
    
    delong_test_result<- roc.test(roc1,roc2,method="delong")
    delong_test_result
    
    
    delong_test_p<-as.numeric(delong_test_result[["p.value"]])
    
    
    delong_test <- rbind(delong_test,data.frame(selected_protein = importance_tune$Feature[i],
                                                Delong_P=delong_test_p))
    
  }
  else{
    
    selected_features_1 <- importance_tune$Feature[1:i]
    
    dat_protein_select_1 <- as.matrix(dat_train[, c("NAFLD",selected_features_1)])
    
    dat_protein_select_train_1 <- lgb.Dataset(data = dat_protein_select_1[, -1], 
                                              label = dat_protein_select_1[, 1])
    
    
    dat_protein_select_test_1<-dat_test[, c("NAFLD",selected_features_1)]
    
    

    lgb_params <- list(objective = 'binary',metric = 'auc',
                       n_estimators=lightgbm_bst_NAFLD$trees,
                       min_data_in_leaf=lightgbm_bst_NAFLD$min_n,
                       max_depth =lightgbm_bst_NAFLD$tree_depth,
                       learning_rate  = lightgbm_bst_NAFLD$learn_rate,
                       min_gain_to_split =lightgbm_bst_NAFLD$loss_reduction)   
    
    
    lgb_model_1 <- lgb.train(
      params = lgb_params,
      data = dat_protein_select_train_1,
      valids = list(valid = lgb.Dataset(data = dat_protein_select_1[, -1], label = dat_protein_select_1[, 1]))
    )
    

    lgb_predicted_prob_tr_1 <- predict(lgb_model_1, as.matrix(dat_protein_select_test_1[,-1]))
    roc1 <- roc(dat_protein_select_test_1$NAFLD, lgb_predicted_prob_tr_1)
    
    selected_features_2 <- importance_tune$Feature[1:i-1]
    
    dat_protein_select_2 <- as.matrix(dat_train[, c("NAFLD",selected_features_2)])
    
    dat_protein_select_train_2 <- lgb.Dataset(data = dat_protein_select_2[, -1], 
                                              label = dat_protein_select_2[, 1])
    
    
    dat_protein_select_test_2<-dat_test[, c("NAFLD",selected_features_2)]
    
    
    lgb_model_2 <- lgb.train(
      params = lgb_params,
      data = dat_protein_select_train_2,
      valids = list(valid = lgb.Dataset(data = dat_protein_select_2[, -1], label = dat_protein_select_2[, 1]))
    )
    

    lgb_predicted_prob_tr_2 <- predict(lgb_model_2, as.matrix(dat_protein_select_test_2[,-1]))
    
    roc2 <- roc(dat_protein_select_test_2$NAFLD, lgb_predicted_prob_tr_2)
    
    delong_test_result<- roc.test(roc1,roc2,method="delong")
    
    
    delong_test_p<-as.numeric(delong_test_result[["p.value"]])
    
    
    delong_test <- rbind(delong_test,data.frame(selected_protein = importance_tune$Feature[i],
                                                Delong_P=delong_test_p))
    
  } 
  
}

write.csv(delong_test,"E:/UKBPPP/Final_result/NAFLD_delong_test.csv")




