#This R script requires pre-installation of three libraries: IRanges in the bioconductor, glmnet and Matrix from CRAN server

#Input
#args[1], path to the folder of tissue-specific enhancer-enhancer and enhancer-gene networks

#args[2], Disease_name: one string provided by users, only used for output file name, can be any name

#args[3], path to the golden standard file
#For GWAS, this file must have three columns as Chr, Position, p_value
#For disease genes, this file should contain a list of known disease genes with their official gene symbols

#args[4], Golden_Standard_flag, a flag indicating whether it is a GWAS file or a gene list file
#For GWAS, Golden_Standard_flag = 1,
#For disease genes, Golden_Standard_flag = 0

#args[5], cutoff_threshold: for GWAS, a cutoff_threshold should be provided, default 0.01

#args[6], path of the output file

#Output
#Output_error_flag
#0: no errors, a file with predicted enhancers is created
#1: error occurs, not enough inputs, program won't run
#2: error occurs, not enough golden standards after mapping, minimum number 20
#3: warning occurs, no good prediction over 0.52

#args <- commandArgs(trailingOnly = TRUE)

#path_tissue_network_folder <- as.character(args[1])
#Disease_name <- as.character(args[2])
#path_golden_standard_file <- as.character(args[3])
#Golden_Standard_flag <- as.numeric(args[4])
#cutoff_threshold <- as.numeric(args[5])
#output_file = args[6]

setwd('~/Desktop/EE_network/FENRIR/')
path_tissue_network_folder <- 'brain'
Disease_name <- 'Disease'
path_golden_standard_file <- 'Example_GWAS_input.txt'
Golden_Standard_flag <- 1
cutoff_threshold <- 0.01
output_file = 'Example_GWAS_output3.txt'

library(IRanges)
library(Matrix)
library(glmnet)
start_time <- Sys.time()

Output_error_flag = 0

# load enhancer regions
Enhancer_region <- read.csv(paste(path_tissue_network_folder,'/Enhancer_regions.txt', sep = ''), header = TRUE, sep = "\t")
Enhancer_region$Chr = as.character(Enhancer_region$Chr)
Num_enhancers = length(Enhancer_region$Enhancer_ID)

# load prior enhancer-gene interactions
EP_network_prior <- read.csv(paste(path_tissue_network_folder,'/Enhancer_gene_network_prior.txt', sep = ''), header = TRUE, sep = "\t")
EP_network_prior_gene_list <- read.csv(paste(path_tissue_network_folder,'/Gene_list.txt', sep = ''), header= TRUE, sep = "\t")
Num_genes = length(as.matrix(EP_network_prior_gene_list))

EP_score=matrix(1, nrow = nrow(EP_network_prior), ncol = 1)
colnames(EP_score) <- 'Score'
EP_network_prior = cbind(EP_network_prior, EP_score)
EP_network_prior_matrix <- Matrix(0, nrow = Num_enhancers, ncol = Num_genes, sparse = TRUE)
EP_network_prior_matrix[cbind(EP_network_prior$Enhancer_ID, EP_network_prior$Gene_ID)] <- EP_network_prior$Score

# load the tissue-specific maximum likely target gene
Enhancer_gene <- read.csv(paste(path_tissue_network_folder,'/Enhancer_gene_network_max_gene.txt', sep = ''), header = TRUE, sep = "\t")
Enhancer_region$Gene <- Enhancer_gene$Gene  
Enhancer_region$Gene = as.character(Enhancer_region$Gene)

# load gene TSS info from hg19
hg19_TSS <- read.csv(paste(path_tissue_network_folder,'/../hg19_TSS.txt', sep=''), header = TRUE, sep = "\t")
index = match(Enhancer_region$Gene, hg19_TSS$Gene_symbol)
Enhancer_region$Gene_chr = hg19_TSS$chr[index]
Enhancer_region$Gene_TSS = hg19_TSS$TSS[index]
Enhancer_region$Gene_strand = hg19_TSS$Strand[index]

# load the tissue-specific enhancer network
EE_network <- read.csv(paste(path_tissue_network_folder,'/Enhancer_Enhancer_network.txt', sep = ''), header = TRUE, sep = "\t")
# transform the three column network to a sparse symmetric matrix
EE_network$Score = log10(EE_network$Score+1)
upper_bound = quantile(EE_network$Score, 0.999)
EE_network$Score = EE_network$Score/upper_bound
EE_network$Score[EE_network$Score>1]=1
EE_network_matrix <- Matrix(0, nrow = Num_enhancers, ncol = Num_enhancers, sparse = TRUE)
EE_network_matrix[cbind(EE_network$Enhancer_ID_1, EE_network$Enhancer_ID_2)] <- EE_network$Score
EE_network_matrix[cbind(EE_network$Enhancer_ID_2, EE_network$Enhancer_ID_1)] <- EE_network$Score

# check if GWAS (1) or diseaase gene file (0) is used
if (Golden_Standard_flag == 1){
  
  # load GWAS file
  GWAS_data <- read.csv(path_golden_standard_file, header = TRUE, sep = "\t")
  GWAS_data$Chr = as.character(GWAS_data$Chr)
  
  GWAS_data$p_value[which(GWAS_data$p_value > cutoff_threshold)] = -1
  GWAS_data$p_value[which(GWAS_data$p_value <= cutoff_threshold & GWAS_data$p_value > 0)] = 1
  
  Enhancer_disease_flag = matrix(0, nrow=Num_enhancers, ncol = 1)
  for (c in unique(Enhancer_region$Chr)){
    enhancer_index=which(Enhancer_region$Chr==c)
    gwas_index=which(GWAS_data$Chr==c)
    for (e in enhancer_index) {
      index=which(GWAS_data$Chr[gwas_index] == Enhancer_region$Chr[e] & 
                    GWAS_data$Position[gwas_index] <=Enhancer_region$Point2[e] & 
                    GWAS_data$Position[gwas_index] >=Enhancer_region$Point1[e])
      if (length(index) > 0){
        Enhancer_disease_flag[e] = max(GWAS_data$p_value[gwas_index[index]])
      }
    }
  }

  #Enhancer_disease_flag_2 = matrix(0, nrow=Num_enhancers, ncol = 1)
  #for (e in c(1:Num_enhancers)) {
  #  index=which(GWAS_data$Chr == Enhancer_region$Chr[e] & 
  #                GWAS_data$Position <=Enhancer_region$Point2[e] & 
  #                GWAS_data$Position >=Enhancer_region$Point1[e])
  #  if (length(index) > 0){
  #    Enhancer_disease_flag_2[e] = max(GWAS_data$p_value[index])
  #  }
  #}


  
  # select positive enhancers
  positive_e_index = which(Enhancer_disease_flag == 1);
  
  # select potential positive target genes of positive enhancers
  positive_g_index_temp = which(colSums(as.matrix(EP_network_prior_matrix[positive_e_index,])) > 0)
  
  # select negative enhancers, non-significant in GWAS and not regulating any potential positive target genes
  negative_e_index = which(Enhancer_disease_flag < 0 & rowSums(as.matrix(EP_network_prior_matrix[,positive_g_index_temp])) == 0)
  
  # check whether get enough positive enhancers, if less than 10, the quite
  if (length(positive_e_index) < 10){
    Output_error_flag = 2
  }
  
}else{
  
  # load user uploaded disease-specific genes
  positive_disease_genes <- read.csv(path_golden_standard_file, header = FALSE)
  
  # load tissue-specific enhancer-gene interactions
  EP_network_max_enhancer <- read.csv(paste(path_tissue_network_folder,'/Enhancer_gene_network_max_enhancer.txt', sep = ''), header = TRUE, sep = "\t")
  
  # select positive enhancers for disease-specific genes
  match_list <- match(EP_network_max_enhancer$Gene, positive_disease_genes[,1])
  positive_e_index = unique(EP_network_max_enhancer$Enhancer_ID[which(match_list > 0)])

  # select potential positive enhancers for disease-specific genes
  match_list <- match(EP_network_prior_gene_list, positive_disease_genes[,1])
  positive_e_index_temp = which(rowSums(as.matrix(EP_network_prior_matrix[,which(match_list > 0)])) > 0)

  # load OMIM genes as background genes and exclude user uploaded disease-sepcific genes
  negative_disease_genes <- read.csv(paste(path_tissue_network_folder,'/../OMIM_genes.txt', sep=''), header = FALSE)
  negative_disease_genes = setdiff(negative_disease_genes, positive_disease_genes)
  
  # select negative enhancers, regulating non-disease-specific genes and not potentially regulating any disease-specific genes
  match_list <- match(EP_network_max_enhancer$Gene, negative_disease_genes[,1])
  negative_e_index = unique(EP_network_max_enhancer$Enhancer_ID[which(match_list > 0)])
  negative_e_index = setdiff(negative_e_index, positive_e_index);
  negative_e_index = setdiff(negative_e_index, positive_e_index_temp);
  
  # check whether we get enough positives
  if (length(positive_e_index) < 10){
    Output_error_flag = 2
  }
}

#Gold standard positive (1) and nagetive (-1) enhancers
Golden_Standard_enhancer_flag = matrix(0, nrow = Num_enhancers, ncol = 1)
Golden_Standard_enhancer_flag[positive_e_index,1] = 1
Golden_Standard_enhancer_flag[negative_e_index,1] = -1

if (Output_error_flag == 2){
  
  #if not enough positives, we won't do training/prediction
  stop("Not enough positive enhancers for training!")  
  
}else{
  
  set.seed(1)
  #round 1 cross validation
  CV_index = c(positive_e_index, negative_e_index)
  CV_labels = matrix(-1, nrow = length(CV_index), ncol = 1)
  CV_labels[c(1:length(positive_e_index))] = 1
  CV_labels = factor(CV_labels)
  random_index = sample.int(length(CV_index))
  CV_labels = CV_labels[random_index];
  CV_index = CV_index[random_index];
  EE_network_matrix_CV = EE_network_matrix[CV_index, ]
  glm_model_1 = cv.glmnet(EE_network_matrix_CV, CV_labels, family = "binomial", alpha = 0.1, type.measure = "auc", nfolds = 5, type.logistic = 'Newton', type.multinomial = 'ungrouped')
  AUC1 = max(glm_model_1$cvm)
  enhancer_probability_1 = predict(glm_model_1, EE_network_matrix, s = 'lambda.min', type = "response")
  probability_sd1 = sd(enhancer_probability_1)
  
  
  #round 2 cross validation
  CV_index = c(positive_e_index, negative_e_index)
  CV_labels = matrix(-1, nrow = length(CV_index), ncol = 1)
  CV_labels[c(1:length(positive_e_index))] = 1
  CV_labels = factor(CV_labels)
  random_index = sample.int(length(CV_index))
  CV_labels = CV_labels[random_index];
  CV_index = CV_index[random_index];
  EE_network_matrix_CV = EE_network_matrix[CV_index, ]
  glm_model_2 = cv.glmnet(EE_network_matrix_CV, CV_labels, family = "binomial", alpha = 0.1, type.measure = "auc", nfolds = 5, type.logistic = 'Newton', type.multinomial = 'ungrouped')
  AUC2 = max(glm_model_2$cvm)
  enhancer_probability_2 = predict(glm_model_2, EE_network_matrix, s = 'lambda.min', type = "response")
  probability_sd2 = sd(enhancer_probability_2)
  
  
  #round 3 cross validation
  CV_index = c(positive_e_index, negative_e_index)
  CV_labels = matrix(-1, nrow = length(CV_index), ncol = 1)
  CV_labels[c(1:length(positive_e_index))] = 1
  CV_labels = factor(CV_labels)
  random_index = sample.int(length(CV_index))
  CV_labels = CV_labels[random_index];
  CV_index = CV_index[random_index];
  EE_network_matrix_CV = EE_network_matrix[CV_index, ]
  glm_model_3 = cv.glmnet(EE_network_matrix_CV, CV_labels, family = "binomial", alpha = 0.1, type.measure = "auc", nfolds = 5, type.logistic = 'Newton', type.multinomial = 'ungrouped')
  AUC3 = max(glm_model_3$cvm)
  enhancer_probability_3 = predict(glm_model_3, EE_network_matrix, s = 'lambda.min', type = "response")
  probability_sd3 = sd(enhancer_probability_3)

  
  #round 4 cross validation
  CV_index = c(positive_e_index, negative_e_index)
  CV_labels = matrix(-1, nrow = length(CV_index), ncol = 1)
  CV_labels[c(1:length(positive_e_index))] = 1
  CV_labels = factor(CV_labels)
  random_index = sample.int(length(CV_index))
  CV_labels = CV_labels[random_index];
  CV_index = CV_index[random_index];
  EE_network_matrix_CV = EE_network_matrix[CV_index, ]
  glm_model_4 = cv.glmnet(EE_network_matrix_CV, CV_labels, family = "binomial", alpha = 0.1, type.measure = "auc", nfolds = 5, type.logistic = 'Newton', type.multinomial = 'ungrouped')
  AUC4 = max(glm_model_4$cvm)
  enhancer_probability_4 = predict(glm_model_4, EE_network_matrix, s = 'lambda.min', type = "response")
  probability_sd4 = sd(enhancer_probability_4)

  
  #round 5 cross validation
  CV_index = c(positive_e_index, negative_e_index)
  CV_labels = matrix(-1, nrow = length(CV_index), ncol = 1)
  CV_labels[c(1:length(positive_e_index))] = 1
  CV_labels = factor(CV_labels)
  random_index = sample.int(length(CV_index))
  CV_labels = CV_labels[random_index];
  CV_index = CV_index[random_index];
  EE_network_matrix_CV = EE_network_matrix[CV_index, ]
  glm_model_5 = cv.glmnet(EE_network_matrix_CV, CV_labels, family = "binomial", alpha = 0.1, type.measure = "auc", nfolds = 5, type.logistic = 'Newton', type.multinomial = 'ungrouped')
  AUC5 = max(glm_model_5$cvm)
  enhancer_probability_5 = predict(glm_model_5, EE_network_matrix, s = 'lambda.min', type = "response")
  probability_sd5 = sd(enhancer_probability_5)
  
  
  best_training=which.max(c(probability_sd1, probability_sd2, probability_sd3, probability_sd4, probability_sd5))
  
  if (best_training==1){
    glm_model=glm_model_1
    enhancer_probability=enhancer_probability_1
  }else if (best_training==2){
    glm_model=glm_model_2
    enhancer_probability=enhancer_probability_2
  }else if (best_training==3){
    glm_model=glm_model_3
    enhancer_probability=enhancer_probability_3
  }else if (best_training==4){
    glm_model=glm_model_4
    enhancer_probability=enhancer_probability_4
  }else{
    glm_model=glm_model_5
    enhancer_probability=enhancer_probability_5
  }

  
  #sort enhancers according to the predicted probability
  sorted_enhancer_probability = sort(enhancer_probability, decreasing = TRUE, index.return = TRUE)
  
  #extract the first 5000 enhancers (about top 10%)
  predicted_enhancers = list('Rank' = c(1:Num_enhancers),
                             'Chr' = Enhancer_region$Chr[sorted_enhancer_probability$ix],
                             'Start' = Enhancer_region$Point1[sorted_enhancer_probability$ix],
                             'End' = Enhancer_region$Point2[sorted_enhancer_probability$ix],
                             'Disease_Probability' = round(1e6*sorted_enhancer_probability$x)/1e6,
                             'Enhancer_Label'= Golden_Standard_enhancer_flag[sorted_enhancer_probability$ix],
                             'Gene_Symbol' = Enhancer_region$Gene[sorted_enhancer_probability$ix],
                             'Gene_Chr' = Enhancer_region$Gene_chr[sorted_enhancer_probability$ix],
                             'Gene_TSS' = Enhancer_region$Gene_TSS[sorted_enhancer_probability$ix],
                             'Gene_Strand' = Enhancer_region$Gene_strand[sorted_enhancer_probability$ix])
  
  # the output file name contains tissue name, disease name, GWAS or gene info
  write.table('*************************************************************************************', file = output_file, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('*** ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(paste('*** ', Disease_name, '-enhancer associations prediction using the ', path_tissue_network_folder, '-specific enhancer network!', sep=''), file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('*** ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  
  if (Golden_Standard_flag == 1){
    write.table(paste('*** Positive enhancers used for training are enhancers containing significant GWAS SNPs with p-value <', cutoff_threshold, sep = ''), file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  }else{
    write.table('*** Positive enhancers used for training are enhancers regulating disease genes', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  }
  
  # Output top predictions with Rank, Chr, Start, End, Disease Probability, Enhancer Label, Gene Symbol, Gene Chr, Gene TSS, Gene Strand
  write.table('*** ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(paste('*** Number of positive enhancers used for training', length(positive_e_index), sep = ': '), file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('*** ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(paste('*** Number of negative enhancers used for training', length(negative_e_index), sep = ': '), file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('*** ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(paste('*** Training AUC', round(1000*max(glm_model$cvm))/1000, sep = ': '), file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('*** ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('*************************************************************************************', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table('  ', file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(paste('Rank','Chr', 'Start', 'End', 'Disease_Probability', 'Enhancer_label', 'Gene_Symbol', 'Gene_Chr', 'Gene_TSS', 'Gene_Strand', sep = '\t'), file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(predicted_enhancers, file = output_file, append = TRUE, col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
  
  #program ends here
}

end_time <- Sys.time()
