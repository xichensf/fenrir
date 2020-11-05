clear
close all
clc
load Enhancer_prior_knowledge
%runing the preprocessing requires a large memery (at least 20 GB)
%% Calculate Pearson Correlation Coeffcient and Z-score for Pairwise Enhancers 
fprintf('Pre-step 1: Calculate Pearson Correlation Coeffcient and Z-score for Pairwise Enhancers\n\n');
EE_cobinding_Z_score=corrcoef(Enhancer_TF_binding_matrix');
EE_cobinding_Z_score(isnan(EE_cobinding_Z_score))=0;
EE_cobinding_Z_score(EE_cobinding_Z_score>0.95)=0.95;
EE_cobinding_Z_score(EE_cobinding_Z_score<-0.95)=-0.95;
EE_cobinding_Z_score=log((1+EE_cobinding_Z_score)./(1-EE_cobinding_Z_score))/2;
EE_cobinding_probability=int8(50*normcdf(EE_cobinding_Z_score, mean(mean(EE_cobinding_Z_score)), 1));
Z_score_vecor=EE_cobinding_Z_score(Enhancer_Enhancer_physical_interactions>0);
EE_cobinding_Z_score_random_1to500_edges_mean = zeros(1, 500);
EE_cobinding_Z_score_random_1to500_edges_std = zeros(1, 500);
Max_edge_num=500;
Num_samples=1000;
for k=1:Max_edge_num
    sample_vector=zeros(1,Num_samples);
    % randomly sample k interactions as null distribution
    for ss=1:Num_samples
        sample_vector(1,ss) = sum(Z_score_vecor(randperm(length(Z_score_vecor), k)))/sqrt(k);
    end   
    EE_cobinding_Z_score_random_1to500_edges_mean(k) = mean(sample_vector);
    EE_cobinding_Z_score_random_1to500_edges_std(k) = std(sample_vector);   
end
EE_cobinding_Z_score=int8(50*EE_cobinding_Z_score);


%% Fit a Poisson dirtibution to ChIA-PET read counts of 3D chromatin interactions
fprintf('Pre-step 2: Fit a Poisson dirtibution to ChIA-PET read counts of 3D chromatin interactions\n\n');
Poisson_lambda=mean(full([Enhancer_Enhancer_ChIAPET_readcount(Enhancer_Enhancer_physical_interactions>0); Enhancer_gene_ChIAPET_readcount(Enhancer_gene_physical_interactions>0)]));
EE_readcount_probability=int8(50*poisscdf(full(Enhancer_Enhancer_ChIAPET_readcount),Poisson_lambda));
EE_readcount_random_1to500_edges_mean = zeros(1, 500);
EE_readcount_random_1to500_edges_std = zeros(1, 500);
EE_read_count_vector=full(Enhancer_Enhancer_ChIAPET_readcount(Enhancer_Enhancer_physical_interactions>0));
for k=1:Max_edge_num
    sample_vector=zeros(1,Num_samples);
    for ss=1:Num_samples
        sample_vector(1,ss) = sum(EE_read_count_vector(randperm(length(EE_read_count_vector), k)))/sqrt(k);
    end   
    EE_readcount_random_1to500_edges_mean(k) = mean(sample_vector);
    EE_readcount_random_1to500_edges_std(k) = std(sample_vector);   
end
EG_readcount_probability=int8(50*poisscdf(full(Enhancer_gene_ChIAPET_readcount),Poisson_lambda));
EG_readcount_random_1to500_edges_mean = zeros(1, 500);
EG_readcount_random_1to500_edges_std = zeros(1, 500);
EG_read_count_vector=full(Enhancer_gene_ChIAPET_readcount(Enhancer_gene_physical_interactions>0));
for k=1:Max_edge_num
    sample_vector=zeros(1,Num_samples);
    for ss=1:Num_samples
        sample_vector(1,ss) = sum(EG_read_count_vector(randperm(length(EG_read_count_vector), k)))/sqrt(k);
    end   
    EG_readcount_random_1to500_edges_mean(k) = mean(sample_vector);
    EG_readcount_random_1to500_edges_std(k) = std(sample_vector);   
end


%% Calculate a joint prior probability of read count and TF binding similirity of pariwise enhancers
fprintf('Pre-step 3: Calculate a joint prior probability of read count and TF binding similirity for each enhancer pair\n\n');
EE_joint_readcount_cobinding_probability=int8(zeros(size(EE_cobinding_Z_score)));
EE_joint_readcount_cobinding_probability(1:10000,:)=int8(sqrt(double(EE_readcount_probability(1:10000,:)).*double(EE_cobinding_probability(1:10000,:))));
EE_joint_readcount_cobinding_probability(10001:20000,:)=int8(sqrt(double(EE_readcount_probability(10001:20000,:)).*double(EE_cobinding_probability(10001:20000,:))));
EE_joint_readcount_cobinding_probability(20001:30000,:)=int8(sqrt(double(EE_readcount_probability(20001:30000,:)).*double(EE_cobinding_probability(20001:30000,:))));
EE_joint_readcount_cobinding_probability(30001:40000,:)=int8(sqrt(double(EE_readcount_probability(30001:40000,:)).*double(EE_cobinding_probability(30001:40000,:))));
EE_joint_readcount_cobinding_probability(40001:end,:)=int8(sqrt(double(EE_readcount_probability(40001:end,:)).*double(EE_cobinding_probability(40001:end,:))));


%% Extract functional relatedness of target genes from GIANT tissue-specific functional gene network 
load GIANT_brain_gene_network.mat
fprintf('Pre-step 4: Extract functional relatedness between target genes from GIANT tissue-specific functional gene network\n\n');
[Target_gene_symbols, bidex, cidex]=intersect(Target_gene_symbols, Gene_symbols, 'stable');
Enhancer_gene_ChIAPET_readcount = Enhancer_gene_ChIAPET_readcount(:,bidex);
Enhancer_gene_physical_interactions = Enhancer_gene_physical_interactions(:,bidex);
GG_functional_relatedness_score=double(GG_functional_relatedness_score(cidex,cidex))/100;
GG_functional_relatedness_score_random_1to500_edges_mean = zeros(1, 500);
GG_functional_relatedness_score_random_1to500_edges_std = zeros(1, 500);
GG_functional_relatedness_score_vector=GG_functional_relatedness_score(Enhancer_gene_physical_interactions'*Enhancer_Enhancer_physical_interactions*Enhancer_gene_physical_interactions>0);
for k=1:Max_edge_num
    sample_vector=zeros(1,Num_samples);
    for ss=1:Num_samples
        sample_vector(1,ss) = sum(GG_functional_relatedness_score_vector(randperm(length(GG_functional_relatedness_score_vector), k)))/sqrt(k);
    end
    GG_functional_relatedness_score_random_1to500_edges_mean(k) = mean(sample_vector);
    GG_functional_relatedness_score_random_1to500_edges_std(k) = std(sample_vector);
end

GG_functional_relatedness_score=int8(50*GG_functional_relatedness_score);

save('FENRIR_input.mat', 'EE_cobinding_probability', 'EE_cobinding_Z_score', 'EE_cobinding_Z_score_random_1to500_edges_mean', 'EE_cobinding_Z_score_random_1to500_edges_std',...
    'EE_joint_readcount_cobinding_probability', 'EE_readcount_probability', 'EE_readcount_random_1to500_edges_mean', 'EE_readcount_random_1to500_edges_std',...
    'EG_readcount_probability', 'EG_readcount_random_1to500_edges_mean', 'EG_readcount_random_1to500_edges_std',...
    'Enhancer_Enhancer_ChIAPET_readcount', 'Enhancer_Enhancer_physical_interactions', 'Enhancer_gene_ChIAPET_readcount', 'Enhancer_gene_physical_interactions',...
    'Enhancer_ID', 'Enhancer_region', 'Enhancer_TF_binding_matrix','Target_gene_symbols', 'TF_symbols', ...
    'GG_functional_relatedness_score', 'GG_functional_relatedness_score_random_1to500_edges_mean', 'GG_functional_relatedness_score_random_1to500_edges_std','-v7.3')

fprintf('FENRIR input data created\n');
