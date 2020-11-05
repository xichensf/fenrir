function [acception_rate, alp]=network_sampling(Enhancer_region, Target_gene_symbols, EE_network, EG_network, EE_read_count, EG_read_count,...
    EE_B_P, EG_P, EE_B_Z_score, GG_functional_score, ...
    EG_ChIA_count_mean, EG_ChIA_count_std, EE_B_Z_score_mean, EE_B_Z_score_std, EE_ChIA_count_mean, EE_ChIA_count_std, GG_functional_score_mean, GG_functional_score_std,...
    Num_enhancers, Num_genes, Max_enhancer_Num, Min_enhancer_num, Max_network_Num, Max_sampling_Num, node_index, tissue_name)
%Metropolis Sampling for FENRIR tissue-specific functional enhancer network inference
%INPUT
% Num_enhancers: the total number of enhancer regions
% Num_genes: the total number of target genes
% EE_network: a binary matrix with 'Num_enhancers' rows and 'Num_enhancers' columns, representing observed physical Enhancer-Enhancer (EE) interactions
% EE_read_count: ChIA-PET read counts matrix corresponds to 'EE_network'
% EE_ChIA_count_mean: a 1 by 500 vector, with each unit as the mean read count for a given number (1 ~ 500) of EE edges  
% EE_ChIA_count_std: a 1 by 500 vector, with each unit as the read count standard deviation for a given number (1 ~ 500) of EE edges 
% EE_B_Z_score: a matrix with 'Num_enhancers' rows and 'Num_genes' columns, each unit as Z-score of TF binding similarity (Pearson Correlation coefficent) for each enhancer pair
% EE_B_Z_score_mean: a 1 by 500 vector, with each unit as the mean Z-score for a given number (1 ~ 500) of EE edges 
% EE_B_Z_score_std: a 1 by 500 vector, with each unit as the Z-score standard deviation for a given number (1 ~ 500) of EE edges 
% EE_B_P: a matrix with 'Num_enhancers' rows and 'Num_enhancers' columns, each unit as a joint probablity of ChIA-PET read count and TF binding similarity for each enhancer pair
% EG_network: a binary matrix with 'Num_enhancers' rows and 'Num_genes' columns, representing observed physical Enhancer-Gene (EG) interactions 
% EG_read_count: ChIA-PET read counts matrix corresponds to 'EG_network'
% EG_ChIA_count_mean: a 1 by 500 vector, with each unit as the mean read count for a given number (1 ~ 500) of EG edges  
% EG_ChIA_count_std: a 1 by 500 vector, with each unit as the read count standard deviation for a given number (1 ~ 500) of EG edges 
% EG_P: a matrix with 'Num_enhancers' rows and 'Num_genes' columns,
% GG_functional_score: a matrix with 'Num_genes' rows and 'Num_genes' columns, each unit as a score or a probablity of functional relatedness for each gene pair
% GG_functional_score_mean: a 1 by 500 vector, with each unit as the mean functional score for a given number (1 ~ 500) of GG edges 
% GG_functional_score_std: a 1 by 500 vector, with each unit as the functional score standard deviation for a given number (1 ~ 500) of EE edges 
% Max_enhancer_Num: the maximum number of enhancers in a sampled enhancer network
% Min_enhancer_num: the minimum number of enhancers in a sampled enhancer network 
% Max_network_Num: the maximum number of EE edges in a sample enhancer network (control network density) 
% Max_sampling_Num: the total number of samppling rounds 
% node_index: specify a node ID to seperate results running on different nodes 
% tissue_name: the name of a human tissue
%OUTPUT
% acception_rate: a 1 by Max_sampling_Num vector, with each unit refelcting given a number of sampling rounds, how many proposals are accepted  
% alp: a 1 by Max_sampling_Num vector, with each unit containing the alpha value calculated for acceptance jude

%% Network initialization

% set the random seed for 
rng('shuffle')

Enhancers = Enhancer_region;
Genes = Target_gene_symbols;

Sampling_network_score = zeros(4, Max_sampling_Num);
alp = zeros(1, Max_sampling_Num);
RATE = zeros(Max_sampling_Num, 4);
proposed_network_score = zeros(Max_sampling_Num, 4);

acception_rate = zeros(1, Max_sampling_Num);
flag = zeros(1, Max_sampling_Num);

EE_network_size = zeros(1, Max_sampling_Num);
EG_network_size = zeros(1, Max_sampling_Num);
GG_network_size = zeros(1, Max_sampling_Num);
EE_MS_statistics = sparse(Num_enhancers, Num_enhancers);
EG_MS_statistics = sparse(Num_enhancers, Num_genes);
GG_MS_statistics = sparse(Num_genes, Num_genes);

% set output timepoint, current very 10% sampling is finished
output_time_point = zeros(10, 1);
for i = 1:9
    output_time_point(i) = round(0.1*i*Max_sampling_Num);   
end
output_time_point(10) = Max_sampling_Num;

% define the upper and lower limit for probabilities
small_data = 0.05;
large_data = 0.95;

%weight for gene functional relatedness
beta_Gene = 2;   

%weight for enhancer co-binding
beta_TF = 2;    

%weight for 3D interaction readcount, weights of both EE and EG types are 2
beta_prior = 1;   

% normalize total weight to 4
beta_Gene = beta_Gene * (4/(beta_Gene + beta_TF + 2*beta_prior));
beta_TF = beta_TF * (4/(beta_Gene + beta_TF+ 2*beta_prior));
beta_prior = beta_prior * (4/(beta_Gene + beta_TF + 2*beta_prior));
total_weight = beta_Gene + beta_TF + 2*beta_prior;

% connectivity degree for all enhancers and select a hub node to start
Enhancer_degree = sum(EE_network);
hub_enhancer_idx = find(Enhancer_degree > 30);
seed_enhancer_idx = hub_enhancer_idx(randperm(length(hub_enhancer_idx), 1));

% select neighbor enhancer nodes
neighbor_nodes = find(EE_network(seed_enhancer_idx,:) == 1);
if length(neighbor_nodes) >= (Max_enhancer_Num - 10)
    neighbor_nodes = neighbor_nodes(randperm(length(neighbor_nodes), Max_enhancer_Num - 10));
end

% construct initial/current EE network
current_E_nodes_idx = union(seed_enhancer_idx, neighbor_nodes);
current_EE_network = sparse(Num_enhancers, Num_enhancers);
current_EE_network(seed_enhancer_idx, neighbor_nodes) = 1;
current_EE_network(neighbor_nodes, seed_enhancer_idx) = 1;
Num_current_EE_interactions = length(find(current_EE_network > 0));
EE_B_Z_score_std(1:end) = 1;

% construct initial EG network
current_EG_network = sparse(Num_enhancers, Num_genes);
for e = 1:length(current_E_nodes_idx)
    gene_index = find(EG_network(current_E_nodes_idx(e), :) > 0);
    current_EG_network(current_E_nodes_idx(e), gene_index(randperm(length(gene_index), 1))) = 1;
end
current_G_nodes_idx = find(sum(current_EG_network) > 0);
Num_current_EG_interactions = length(find(current_EG_network > 0));

% construct initial GG network
% two target genes are connected if their upstream enhancers have a direct interaction
current_GG_network = current_EG_network'*current_EE_network*current_EG_network;
current_GG_network(current_GG_network > 0) = 1;
Num_current_GG_interactions = length(find(current_GG_network > 0));
GG_functional_score_std(1:end) = 1;
GG_P = GG_functional_score;

% calcluate network probability P(X) for the initial network state
current_network_score_EE_prior = (full(sum(EE_read_count(current_EE_network > 0))) / sqrt(Num_current_EE_interactions) - EE_ChIA_count_mean(Num_current_EE_interactions))...
                                  /EE_ChIA_count_std(Num_current_EE_interactions);
current_network_score_EE_B = (sum(EE_B_Z_score(current_EE_network > 0))*0.02 / sqrt(Num_current_EE_interactions) - EE_B_Z_score_mean(Num_current_EE_interactions))...
                              /EE_B_Z_score_std(Num_current_EE_interactions);
current_network_score_EG_prior = (full(sum(EG_read_count(current_EG_network > 0))) / sqrt(Num_current_EG_interactions) - EG_ChIA_count_mean(Num_current_EG_interactions))...
                                  /EG_ChIA_count_std(Num_current_EG_interactions);
current_network_score_GG_function = (sum(GG_functional_score(current_GG_network > 0))*0.02 / sqrt(Num_current_GG_interactions) - GG_functional_score_mean(Num_current_GG_interactions))...
                                     /GG_functional_score_std(Num_current_GG_interactions);
current_network_score = [current_network_score_EE_prior, current_network_score_EE_B, current_network_score_EG_prior, current_network_score_GG_function];


%% ********** Metropolis sampling starts! *************************
for it = 1:Max_sampling_Num
    % generate a random number x to control add or delete actions
    % x>0.5 add; x<0.5 delete
    x = rand;    
    
    % find all EE interactions whose deletion will not break the current EE network
    [EE_xx, EE_yy] = find(current_EE_network>0);
    cutting_edge_flag = -1*ones(length(EE_xx), 1);
    for k = 1:length(EE_xx)
        if EE_xx(k) > EE_yy(k)
            test_network = current_EE_network;
            test_network(EE_xx(k), EE_yy(k)) = 0;           
            test_network(EE_yy(k), EE_xx(k)) = 0;
            current_eidx = find(sum(test_network) > 0);
            A = full(test_network(current_eidx, current_eidx));
            D = diag(sum(A));
            L = D - A;
            eigen_value = eig(L);
            eigen_value = sort(eigen_value, 'ascend');
            cutting_edge_flag(k) = eigen_value(2);
        end
    end  
    
    % if this edge is deleted, the network is still connected, second largest egen value >0
    cutting_edge_index = find(cutting_edge_flag > 0.001);

    % add action
    if ((x >= 0.5 && length(current_E_nodes_idx) < Max_enhancer_Num) && (Num_current_EE_interactions < Max_network_Num))...
            || (length(current_E_nodes_idx) <= Min_enhancer_num)...
            || (length(cutting_edge_index) <= 5)
        
        flag(it) = 1;
        
        % propose a network state X' by adding an EE edge and an EG edge
        % calculate the conditional probability P(X'|X)
        [proposed_EE_network, proposed_EG_network, proposed_GG_network,...
         seed_enhancer_idx, added_Enhancer_node_idx, seed_gene_idx, added_Gene_node_idx,...
         adding_enhancer_prior, adding_enhancer_probability, adding_gene_prior, adding_gene_probability] = ...
                add_jump_function(current_EE_network, current_EG_network, current_GG_network,EE_network, EG_network, EE_B_P, EG_P, GG_P,...
                                  small_data, large_data, Num_enhancers, Max_enhancer_Num);
        
        % reversable jump 
        % calculate the conditional probability P(X|X')
        [reverse_delete_prior, reverse_delete_probability] =...
                reverse_delete_jump_function(proposed_EE_network, proposed_EG_network,seed_enhancer_idx, added_Enhancer_node_idx, seed_gene_idx, added_Gene_node_idx,...
                                             EE_B_P, EE_B_P, EG_P, GG_P,small_data, large_data);
        
        % calculate P(X|X')/P(X|X')
        conditional_prior_ratio = (reverse_delete_probability/reverse_delete_prior)...
                                   /sqrt((adding_enhancer_probability/adding_enhancer_prior) * (adding_gene_probability/adding_gene_prior));
    
    else
        % delete action
        flag(it) = -1;
        
        % propose a network state X' by deleting an EE edge and their target genes in EG
        % calculate the conditional probability P(X'|X)
        [proposed_EE_network, proposed_EG_network, proposed_GG_network,...
         deleted_enhancer_hub_node, deleted_enhancer_leaf_node, deleted_hub_gene, deleted_leaf_gene,...
         delete_edge_prior, delete_edge_probability] = ...
                delete_jump_function(current_EE_network, current_EG_network, current_GG_network, EE_B_P, EE_B_P, EG_P, GG_P,...
                                     small_data, large_data, EE_xx, EE_yy, cutting_edge_index);
        
        % reverse jump
        % calculate the conditional probability P(X|X')
        [reverse_adding_enhancer_prior, reverse_adding_enhancer_probability, reverse_adding_gene_prior, reverse_adding_gene_probability] = ...
                reverse_add_jump_function(proposed_EE_network, proposed_EG_network, proposed_GG_network, EE_network, EG_network, ...
                                          deleted_enhancer_hub_node, deleted_enhancer_leaf_node, deleted_hub_gene, deleted_leaf_gene, EE_B_P, EG_P, GG_P,...
                                          small_data, large_data, Max_enhancer_Num, Num_enhancers);
        
        % calculate P(X|X')/P(X|X')
        conditional_prior_ratio=sqrt((reverse_adding_enhancer_probability/reverse_adding_enhancer_prior)*(reverse_adding_gene_probability/reverse_adding_gene_prior))...
                                    /(delete_edge_probability/delete_edge_prior);

    end
    
    Num_proposed_EE_interactions = sum(proposed_EE_network(proposed_EE_network > 0));
    Num_proposed_EG_interactions = sum(proposed_EG_network(proposed_EG_network > 0));
    Num_proposed_GG_interactions = sum(proposed_GG_network(proposed_GG_network > 0));
    
    % calculate network probability P(X') for the proposed network state
    proposed_network_score_EE_prior = (full(sum(EE_read_count(proposed_EE_network > 0)))/sqrt(Num_proposed_EE_interactions)/sqrt(2) - EE_ChIA_count_mean(Num_proposed_EE_interactions/2))...
                                       /EE_ChIA_count_std(Num_proposed_EE_interactions);
    proposed_network_score_EE_B = (sum(double(EE_B_Z_score(proposed_EE_network > 0))*0.02)/sqrt(Num_proposed_EE_interactions)/sqrt(2) - EE_B_Z_score_mean(Num_proposed_EE_interactions/2))...
                                   /EE_B_Z_score_std(Num_proposed_EE_interactions);
    proposed_network_score_EG_prior = (full(sum(EG_read_count(proposed_EG_network > 0)))/sqrt(Num_proposed_EG_interactions) - EG_ChIA_count_mean(Num_proposed_EG_interactions))...
                                       /EG_ChIA_count_std(Num_proposed_EG_interactions);
    proposed_network_score_GG_function = (sum(double(GG_functional_score(proposed_GG_network > 0))*0.02)/sqrt(Num_proposed_GG_interactions)/sqrt(2) - GG_functional_score_mean(Num_proposed_GG_interactions/2))...
                                          /GG_functional_score_std(Num_proposed_GG_interactions);
    proposed_network_score(it,:) = [proposed_network_score_EE_prior, proposed_network_score_EE_B, proposed_network_score_EG_prior, proposed_network_score_GG_function];
   
    % calculate P(X')/P(X)
    RATE(it,:) = normcdf(proposed_network_score(it,:))./normcdf(current_network_score);
    RATE(it,RATE(it,:) > 5) = 5;
    RATE(it,RATE(it,:) < 0.2) = 0.2;
    
    % calulate alpha = P(X')P(X|X')/(P(X)P(X'|X))
    alp(it) = (RATE(it,1)^beta_prior) * (RATE(it,2)^beta_TF) * (RATE(it,3)^beta_prior) * (RATE(it,4)^beta_Gene) * (conditional_prior_ratio^total_weight);
    
    % judge if the proposed network state X' can be accepted
    temp = rand;
    if (alp(it) > temp || isnan(alp(it)))
        current_network_score = proposed_network_score(it,:);
        current_EE_network = proposed_EE_network;
        current_EG_network = proposed_EG_network;
        current_GG_network = proposed_GG_network; 
        current_E_nodes_idx = find(sum(current_EE_network) > 0);
        current_G_nodes_idx = find(sum(current_EG_network) > 0);
        Num_current_EE_interactions = Num_proposed_EE_interactions;
        Num_current_EG_interactions = Num_proposed_EG_interactions;
        Num_current_GG_interactions = Num_proposed_GG_interactions;
    end
    
    % record network state samples
    EE_MS_statistics = EE_MS_statistics + current_EE_network;
    EG_MS_statistics = EG_MS_statistics + current_EG_network;
    GG_MS_statistics = GG_MS_statistics + current_GG_network; 
    
    Sampling_network_score(:,it) = current_network_score';
    
    if it == 1
        acception_rate(it) = min(alp(it),1);
    else
        acception_rate(it) = (acception_rate(it-1)*(it-1)+min(alp(it),1))/it;
    end
    
    EE_network_size(it) = Num_current_EE_interactions;
    EG_network_size(it) = Num_current_EG_interactions;
    GG_network_size(it) = Num_current_GG_interactions; 
    
    % when finishing every 10% sampling work, output temp results
    if ~isempty(find(output_time_point == it) > 0)
        fprintf([num2str(it), ' rounds of sampling done!\n'])
        save([tissue_name, '_node_', num2str(node_index), '_', num2str(it), '.mat'], 'it', 'Enhancers', 'Genes', 'EE_MS_statistics', 'EG_MS_statistics')
        EE_MS_statistics = sparse(Num_enhancers, Num_enhancers);
        EG_MS_statistics = sparse(Num_enhancers, Num_genes);
        GG_MS_statistics = sparse(Num_genes, Num_genes);
    end    
end

% after finishing all sampling work, output final results
EE_MS_final_statistics = sparse(Num_enhancers, Num_enhancers);
EG_MS_final_statistics = sparse(Num_enhancers, Num_genes);
for i=1:10
    load([tissue_name, '_node_', num2str(node_index), '_', num2str(output_time_point(i)), '.mat'], 'it', 'EE_MS_statistics', 'EG_MS_statistics')
    EE_MS_final_statistics = EE_MS_final_statistics + EE_MS_statistics;
    EG_MS_final_statistics = EG_MS_final_statistics + EG_MS_statistics;
end
save([tissue_name, '_node_', num2str(node_index), '_final_network.mat'], 'Enhancers', 'Genes', 'EE_MS_statistics', 'EG_MS_statistics')

