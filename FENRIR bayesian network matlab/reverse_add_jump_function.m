function [reverse_adding_enhancer_prior, reverse_adding_enhancer_probability, reverse_adding_gene_prior, reverse_adding_gene_probability]=...
    reverse_add_jump_function(proposed_EE_network, proposed_EG_network, proposed_GG_network, ...
                                EE_network, EG_network, ...
                                deleted_enhancer_hub_node, deleted_enhancer_leaf_node, deleted_hub_gene, deleted_leaf_gene,...
                                EE_B_P, EG_P, GG_P,...
                                small_data, large_data, Max_enhancer_Num, Num_enhancers)    
    current_neighbor_enhancers=union(deleted_enhancer_hub_node, find(proposed_EE_network(deleted_enhancer_hub_node,:)>0));% find its neighboring enahncers in current network  
    % select Max_network_Num neighbor enhancers with highest EE prior probability    
    [aa, idx]=sort(double(EE_B_P(deleted_enhancer_hub_node, :))+0.5*rand(1,Num_enhancers), 'descend');  
    reverse_neighbor_enhancers=idx(1:2*Max_enhancer_Num);
    
    reverse_neighbor_enhancers=setdiff(reverse_neighbor_enhancers, current_neighbor_enhancers);% select new enhancers
    reverse_adding_enhancer_prior=1/length(reverse_neighbor_enhancers);% prior pdf for all neighbors
    
    reverse_neighbor_enhancers_weight=double(EE_B_P(deleted_enhancer_hub_node, reverse_neighbor_enhancers))*0.02;
    reverse_neighbor_enhancers_weight(reverse_neighbor_enhancers_weight>large_data)=large_data;
    reverse_neighbor_enhancers_weight(reverse_neighbor_enhancers_weight<small_data)=small_data;

    reverse_adding_enhancer_probability=double(EE_B_P(deleted_enhancer_hub_node, deleted_enhancer_leaf_node))*0.02;
    reverse_adding_enhancer_probability(reverse_adding_enhancer_probability>large_data)=large_data;
    reverse_adding_enhancer_probability(reverse_adding_enhancer_probability<small_data)=small_data;    
    reverse_adding_enhancer_probability=reverse_adding_enhancer_probability/sum(reverse_neighbor_enhancers_weight);
        
    % select neighbor promoters with at least one ChIA-PET connection
    reverse_neighbor_genes=find(EG_network(deleted_enhancer_leaf_node,:)>0);% select target genes for the new neighbor enhancer
    reverse_adding_gene_prior=1/length(reverse_neighbor_genes);% prior for all candidate new target genes
    
    if length(reverse_neighbor_genes)==1
        reverse_adding_gene_probability=1;
    elseif ~isempty(find(proposed_EE_network(deleted_enhancer_leaf_node,:)>0)) % if enhancer is already within network, we do not add genes
        reverse_adding_gene_prior=1;
        reverse_adding_gene_probability=1;
    else % only for new enhancer node
        reverse_neighbor_genes_weight=0.02*sqrt(double(EG_P(deleted_enhancer_leaf_node, reverse_neighbor_genes)).*double(GG_P(deleted_hub_gene, reverse_neighbor_genes)));    % sample a node according to its prior probability and co-binding similarity       
        reverse_neighbor_genes_weight(reverse_neighbor_genes_weight>large_data)=large_data;  
        reverse_neighbor_genes_weight(reverse_neighbor_genes_weight<small_data)=small_data;        
        reverse_adding_gene_probability=0.02*sqrt(double(EG_P(deleted_enhancer_leaf_node, deleted_leaf_gene))*double(GG_P(deleted_hub_gene, deleted_leaf_gene)));% sample a node according to its prior probability and co-binding similarity       
        
        if reverse_adding_gene_probability>large_data
            reverse_adding_gene_probability=large_data;
        end
        if reverse_adding_gene_probability<small_data
            reverse_adding_gene_probability=small_data;
        end
        reverse_adding_gene_probability=reverse_adding_gene_probability/sum(reverse_neighbor_genes_weight);
    end
