function [proposed_EE_network, proposed_EG_network, proposed_GG_network, ...
        deleted_enhancer_hub_node, deleted_enhancer_leaf_node, deleted_hub_gene, deleted_leaf_gene,...
        delete_edge_prior, delete_edge_probability]=...
        delete_jump_function(current_EE_network, current_EG_network, current_GG_network,...
                        EE_P, B_P, EG_P, GG_P,...
                        small_data, large_data, EE_xx, EE_yy, cutting_edge_index)
                    
    EE_xx=EE_xx(cutting_edge_index);
    EE_yy=EE_yy(cutting_edge_index);
    
    Num_candidate_edges=length(EE_xx);
    delete_edge_prior=1/Num_candidate_edges;
    
    if Num_candidate_edges==1
        delete_edge_probability=1;
        if length(find(current_EE_network(EE_xx(1),:)>0))>=length(find(current_EE_network(EE_yy(1),:)>0))
            deleted_enhancer_hub_node=EE_xx(1);
            deleted_enhancer_leaf_node=EE_yy(1);
            deleted_hub_gene=find(current_EG_network(deleted_enhancer_hub_node,:)>0);
            deleted_leaf_gene=find(current_EG_network(deleted_enhancer_leaf_node,:)>0);
        else
            deleted_enhancer_hub_node=EE_yy(1);
            deleted_enhancer_leaf_node=EE_xx(1);
            deleted_hub_gene=find(current_EG_network(deleted_enhancer_hub_node,:)>0);
            deleted_leaf_gene=find(current_EG_network(deleted_enhancer_leaf_node,:)>0);
        end
    else       
        candidate_enhancer_edge_weights=zeros(1,Num_candidate_edges);
        for e=1:Num_candidate_edges
            e_hub_idx=EE_xx(e);
            e_leaf_idx=EE_yy(e);
            g_hub_idx=find(current_EG_network(e_hub_idx,:)>0);
            g_leaf_idx=find(current_EG_network(e_leaf_idx,:)>0);
            
            P1=double(EE_P(e_hub_idx,e_leaf_idx))*0.02;
            P2=double(B_P(e_hub_idx,e_leaf_idx))*0.02;       
            P3_1=double(EG_P(e_hub_idx,g_hub_idx(1)))*0.02;
            P3_2=double(EG_P(e_leaf_idx,g_leaf_idx(1)))*0.02;
            P3=sqrt(P3_1*P3_2);
            P4=double(GG_P(g_leaf_idx(1), g_hub_idx(1)))*0.02;         
            candidate_enhancer_edge_weights(1,e)=1-(P1*P2*P3*P4)^0.25;     
        end

        candidate_enhancer_edge_weights(candidate_enhancer_edge_weights>large_data)=large_data;
        candidate_enhancer_edge_weights(candidate_enhancer_edge_weights<small_data)=small_data;
       
        deleted_edge_index=randsample(Num_candidate_edges, 1, true, candidate_enhancer_edge_weights);
        delete_edge_probability=candidate_enhancer_edge_weights(deleted_edge_index)/sum(candidate_enhancer_edge_weights);
        
        % identify hub enhancer and see enhancer for reverse jump
        if length(find(current_EE_network(EE_xx(deleted_edge_index),:)>0))>=length(find(current_EE_network(EE_yy(deleted_edge_index),:)>0))
            deleted_enhancer_hub_node=EE_xx(deleted_edge_index);
            deleted_enhancer_leaf_node=EE_yy(deleted_edge_index);
            deleted_hub_gene=find(current_EG_network(deleted_enhancer_hub_node,:)>0);
            deleted_leaf_gene=find(current_EG_network(deleted_enhancer_leaf_node,:)>0);
        else
            deleted_enhancer_hub_node=EE_yy(deleted_edge_index);
            deleted_enhancer_leaf_node=EE_xx(deleted_edge_index);
            deleted_hub_gene=find(current_EG_network(deleted_enhancer_hub_node,:)>0);
            deleted_leaf_gene=find(current_EG_network(deleted_enhancer_leaf_node,:)>0);
        end

    end
    % delete this EE interaction   
    proposed_EE_network=current_EE_network;
    proposed_EE_network(deleted_enhancer_leaf_node, deleted_enhancer_hub_node)=0;
    proposed_EE_network(deleted_enhancer_hub_node, deleted_enhancer_leaf_node)=0;
    
    % delete genes if the upstream enhancer is not connected with any other enhancer
    proposed_EG_network=current_EG_network;
    if sum(proposed_EE_network(deleted_enhancer_hub_node,:))==0
        proposed_EG_network(deleted_enhancer_hub_node, :)=0;
    end
    if sum(proposed_EE_network(deleted_enhancer_leaf_node,:))==0
        proposed_EG_network(deleted_enhancer_leaf_node, :)=0;
    end
    proposed_GG_network=proposed_EG_network'*proposed_EE_network*proposed_EG_network;
