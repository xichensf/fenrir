function [reverse_delete_prior, reverse_delete_probability] =...
    reverse_delete_jump_function(proposed_EE_network, proposed_EG_network,...
                                seed_enhancer_idx, added_Enhancer_node_idx, seed_gene_idx, added_Gene_node_idx,...
                                EE_P, B_P, EG_P, GG_P,...
                                small_data, large_data) 
    
    [reverse_EE_xx,reverse_EE_yy]=find(proposed_EE_network>0);% selected edges can be deleted from the proposed network
    reverse_cutting_edge_flag=-1*ones(length(reverse_EE_xx), 1);
    for k=1:length(reverse_EE_xx)
        if reverse_EE_xx(k)>reverse_EE_yy(k)
            test_network=proposed_EE_network;
            test_network(reverse_EE_xx(k), reverse_EE_yy(k))=0;
            test_network(reverse_EE_yy(k), reverse_EE_xx(k))=0;
            current_eidx=find(sum(test_network)>0);
            A=full(test_network(current_eidx, current_eidx));%adjacent graph
            D=diag(sum(A));%degree matrix
            L=D-A;
            eigen_value=eig(L);
            eigen_value=sort(eigen_value, 'ascend');
            reverse_cutting_edge_flag(k)=eigen_value(2);
        end
    end  
    
    reverse_cutting_edge_index=find(reverse_cutting_edge_flag>0.001);% if this edge is deleted, the network is still connected, second largest egen value >0
    reverse_Num_candidate_edges=length(reverse_cutting_edge_index);
    reverse_delete_prior=1/reverse_Num_candidate_edges;% prior pdf for reverse delete
    
    reverse_candidate_enhancer_edge_weights=zeros(1, reverse_Num_candidate_edges);
    for e=1:reverse_Num_candidate_edges
        e_hub_idx=reverse_EE_xx(e);
        e_leaf_idx=reverse_EE_yy(e);
        g_hub_idx=find(proposed_EG_network(e_hub_idx,:)>0);            
        g_leaf_idx=find(proposed_EG_network(e_leaf_idx,:)>0);
        
        P1=double(EE_P(e_hub_idx,e_leaf_idx))*0.02;
        P2=double(B_P(e_hub_idx,e_leaf_idx))*0.02;       
        P3_1=double(EG_P(e_hub_idx,g_hub_idx(1)))*0.02;
        P3_2=double(EG_P(e_leaf_idx,g_leaf_idx(1)))*0.02;
        P3=sqrt(P3_1*P3_2);
        P4=double(GG_P(g_leaf_idx(1), g_hub_idx(1)))*0.02;         
        reverse_candidate_enhancer_edge_weights(1,e)=1-(P1*P2*P3*P4)^0.25;
    end
    
    reverse_candidate_enhancer_edge_weights(reverse_candidate_enhancer_edge_weights>large_data)=large_data;
    reverse_candidate_enhancer_edge_weights(reverse_candidate_enhancer_edge_weights<small_data)=small_data;
    
    P1=double(EE_P(seed_enhancer_idx,added_Enhancer_node_idx))*0.02;
	P2=double(B_P(seed_enhancer_idx,added_Enhancer_node_idx))*0.02;       
	P3_1=double(EG_P(seed_enhancer_idx,seed_gene_idx))*0.02;
	P3_2=double(EG_P(added_Enhancer_node_idx,added_Gene_node_idx))*0.02;
	P3=sqrt(P3_1*P3_2);
	P4=double(GG_P(seed_gene_idx, added_Gene_node_idx))*0.02;         
    reverse_delete_probability=1-(P1*P2*P3*P4)^0.25;
    
    if reverse_delete_probability>large_data
        reverse_delete_probability=large_data;
    end
    
    if reverse_delete_probability<small_data
        reverse_delete_probability=small_data;
    end
    reverse_delete_probability=reverse_delete_probability/sum(reverse_candidate_enhancer_edge_weights);% candidational pdf for reverse jump
