clear
close all
clc
% this input file can be generated using the Prior_data_preprocessing.m file   
load FENRIR_input_demo.mat

tissue_name = 'Brain';
% rounds of Metropolis-Hastlings Sampling (i.e. 10000000)
Max_sampling_Num = 1000;   

% upper limit of enhancers in a sampled subnetwork
Max_enhancer_Num = 50; 

% lower limit of enhancers in a sampled subnetwork
Min_enhancer_num = 10;  

% upper limit of EE interactions in a sampled subnetwork
Max_network_Num = 500;    

% the number of enhancers
Num_enhancers = length(Enhancer_ID);

% the number of target genes
Num_genes = length(Target_gene_symbols);

% if this main function was called on different CPU nodes, this node index
% should be unique to avoid output filename conflict
node_index = 1;

% This function will autimatically output sampling statitics when finishing every 10% of Max_sampling_Num and then a final network 
network_sampling(Enhancer_region, Target_gene_symbols, Enhancer_Enhancer_physical_interactions, Enhancer_gene_physical_interactions,...
    Enhancer_Enhancer_ChIAPET_readcount, Enhancer_gene_ChIAPET_readcount,...
    EE_joint_readcount_cobinding_probability, EG_readcount_probability, EE_cobinding_Z_score, GG_functional_relatedness_score,...
    EG_readcount_random_1to500_edges_mean, EG_readcount_random_1to500_edges_std,...
    EE_cobinding_Z_score_random_1to500_edges_mean, EE_cobinding_Z_score_random_1to500_edges_std,...
    EE_readcount_random_1to500_edges_mean, EE_readcount_random_1to500_edges_std,...
    GG_functional_relatedness_score_random_1to500_edges_mean, GG_functional_relatedness_score_random_1to500_edges_std,...
    Num_enhancers, Num_genes, Max_enhancer_Num, Min_enhancer_num, Max_network_Num, Max_sampling_Num, node_index, tissue_name);

