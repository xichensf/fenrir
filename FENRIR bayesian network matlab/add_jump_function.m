function [proposed_EE_network, proposed_EG_network, proposed_GG_network, ...
    seed_enhancer_idx, added_Enhancer_node_idx, seed_gene_idx, added_Gene_node_idx,...
    adding_enhancer_prior, adding_enhancer_probability, adding_gene_prior, adding_gene_probability] =...
    add_jump_function(current_EE_network, current_EG_network, current_GG_network,...
                        EE_network, EG_network, ...
                        EE_B_P, EG_P, GG_P,...
                        small_data, large_data, Num_enhancers, Max_enhancer_Num) 

    current_E_nodes_idx=find(sum(current_EE_network)>0);
    seed_enhancer_idx=current_E_nodes_idx(randperm(length(current_E_nodes_idx), 1));% randomly select a seed enhancer
    current_neighbor_enhancers=union(seed_enhancer_idx, find(current_EE_network(seed_enhancer_idx,:)>0));% find its neighboring enahncers in current network
    
    % select Max_network_Num neighbor enhancers with highest EE prior probability and Co-binding similarity 
    [aa, idx]=sort(double(EE_B_P(seed_enhancer_idx, :))+0.5*rand(1,Num_enhancers), 'descend');  
    candidate_neighbor_enhancers=idx(1:2*Max_enhancer_Num);
    candidate_neighbor_enhancers=setdiff(candidate_neighbor_enhancers, current_neighbor_enhancers);% at least one new node
    adding_enhancer_prior=1/length(candidate_neighbor_enhancers);
    
    candidate_neighbor_enhancers_weight=double(EE_B_P(seed_enhancer_idx, candidate_neighbor_enhancers))*0.02;
    candidate_neighbor_enhancers_weight(candidate_neighbor_enhancers_weight>large_data)=large_data;
    candidate_neighbor_enhancers_weight(candidate_neighbor_enhancers_weight<small_data)=small_data;
%     candidate_neighbor_enhancers_weight=candidate_neighbor_enhancers_weight./(1-candidate_neighbor_enhancers_weight);
    %sample a node according to its prior probability and co-binding similarity
    E_idx= randsample(length(candidate_neighbor_enhancers),1,true,candidate_neighbor_enhancers_weight);
    added_Enhancer_node_idx=candidate_neighbor_enhancers(E_idx);
    adding_enhancer_probability=candidate_neighbor_enhancers_weight(E_idx)/sum(candidate_neighbor_enhancers_weight);
    
    proposed_EE_network=current_EE_network;
    proposed_EE_network(seed_enhancer_idx, added_Enhancer_node_idx)=1;
    proposed_EE_network(added_Enhancer_node_idx, seed_enhancer_idx)=1;
    
    % Step 2, add a target gene
    seed_gene_idx=find(current_EG_network(seed_enhancer_idx, :)==1);% select target gene of the seed enhancer
    
    % select neighbor promoters with at least one ChIA-PET connection
    candidate_neighbor_genes=find(EG_network(added_Enhancer_node_idx,:)==1);% select target genes for the new neighbor enhancer
    adding_gene_prior=1/length(candidate_neighbor_genes);% prior for all candidate new target genes
    
    if length(candidate_neighbor_genes)==1
        added_Gene_node_idx=candidate_neighbor_genes(1);
        adding_gene_probability=1;
    elseif ~isempty(intersect(added_Enhancer_node_idx, current_E_nodes_idx)) % if enhancer is already within network, we do not add genes
        added_Gene_node_idx=find(current_EG_network(added_Enhancer_node_idx, :)>0);
        adding_gene_prior=1;
        adding_gene_probability=1;
    else % only for new enhancer node
        candidate_neighbor_genes_weight=0.02*sqrt(double(EG_P(added_Enhancer_node_idx, candidate_neighbor_genes)).*double(GG_P(seed_gene_idx, candidate_neighbor_genes)));    % sample a node according to its prior probability and co-binding similarity       
        candidate_neighbor_genes_weight(candidate_neighbor_genes_weight>large_data)=large_data;  
        candidate_neighbor_genes_weight(candidate_neighbor_genes_weight<small_data)=small_data; 
        gene_idx = randsample(length(candidate_neighbor_genes),1,true,candidate_neighbor_genes_weight);
        added_Gene_node_idx=candidate_neighbor_genes(gene_idx);
        adding_gene_probability=candidate_neighbor_genes_weight(gene_idx)/sum(candidate_neighbor_genes_weight);% canditional pdf for the new target gene
    end
    proposed_EG_network=current_EG_network;
    proposed_EG_network(added_Enhancer_node_idx, added_Gene_node_idx)=1;
    
    proposed_GG_network=proposed_EG_network'*proposed_EE_network*proposed_EG_network;
