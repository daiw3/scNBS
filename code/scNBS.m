function [mag_effect_pos_final, mag_effect_neg_final, Z_score, edges_pos_final, edges_neg_final] = scNBS(mat,y, colnames, node_to_network, r_method)
    %% split data into 2 sets
    n = size(mat, 3);
    sample1=randsample(n,floor(n/2));
    sample2=setdiff(1:n, sample1);
    sample2=sample2';
    n_networks=histc(node_to_network.network,unique(node_to_network.network));
    
    
    %% Get Ranking for all Blocks/Networks
    %fprintf("Ranking Method: %s.\n", r_method);
    Z_score = scNBS_ranking(mat,y,sample1, colnames, node_to_network, r_method);
    %disp("Ranking done.");
    
    mag_effect_pos_final = cell(length(n_networks), length(n_networks));
    mag_effect_neg_final = cell(length(n_networks), length(n_networks));
    edges_pos_final = cell(length(n_networks), length(n_networks));
    edges_neg_final = cell(length(n_networks), length(n_networks));
    for net1 = 1:length(n_networks)
        for net2 = net1:length(n_networks)
            node_idx1= node_to_network.node(find(node_to_network.network == net1));
            node_idx2= node_to_network.node(find(node_to_network.network == net2));
            X_subnetwork = mat(node_idx1, node_idx2,:);
            
            %% Cut-off Selection
            [cut_pos, cut_neg, pos_edges, neg_edges] = scNBS_cutoff(X_subnetwork,y,sample1, Z_score{net1,net2});
            
             %% Network-level Marginal Testing
             
             %%% Positive
             if ~isempty(cut_pos)
                 X_subnetwork_pos = X_subnetwork(pos_edges(:,2), pos_edges(:,3), sample2);
                 X_subnetwork_pos_vec = zeros([length(sample2),1]);
                 for i = 1:length(sample2)
                     X_subnetwork_pos_vec(i,:) = mean(diag(X_subnetwork_pos(:,:,i)));
                 end

                 m_pos = fitlm(X_subnetwork_pos_vec,y(sample2,:));
                 mag_effect_pos_final{net1,net2} = [table2array(m_pos.Coefficients(2,:))];
                 edges_pos_final{net1, net2} = pos_edges(:,2:3);
             else
                 mag_effect_pos_final{net1, net2} = NaN(1, 4);
             end  
             %%% Negative
             if ~isempty(cut_neg)
                 X_subnetwork_neg = X_subnetwork(neg_edges(:,2), neg_edges(:,3), sample2);
                 X_subnetwork_neg_vec = zeros([length(sample2),1]);
                 for i = 1:length(sample2)
                     X_subnetwork_neg_vec(i,:) = mean(diag(X_subnetwork_neg(:,:,i)));
                 end

                 m_neg = fitlm(X_subnetwork_neg_vec,y(sample2,:));
                 mag_effect_neg_final{net1,net2} = [table2array(m_neg.Coefficients(2,:))];
                 edges_neg_final{net1, net2} = neg_edges(:,2:3);
             else
                 mag_effect_neg_final{net1, net2} = NaN(1, 4);
             end 
             % fprintf('Network Done: %d-%d.\n', net1, net2);
        end
    end
end