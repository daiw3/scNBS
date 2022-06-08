function [cut_best_pos, cut_best_neg, pos_edges, neg_edges] = scNBS_cutoff(mat_net,y,sample1, Z_mat)
    n = size(mat_net, 3);

    %% convert 3D array to 2D matrix
    if size(mat_net,1) == size(mat_net,2)
        X_subnetwork_vec = zeros([n,size(mat_net,1)*(size(mat_net,1)-1)/2]);
        for j=1:n
            tmp = mat_net(:,:,j);
            X_subnetwork_vec(j,:)=tmp(tril(true(size(tmp)),-1));
        end
    else
        X_subnetwork_vec = zeros([n,size(mat_net,1)*size(mat_net,2)]);
        for j=1:n
            tmp=mat_net(:,:,j);
            X_subnetwork_vec(j,:)=tmp(:);
        end 
    end
    X_sub = X_subnetwork_vec(sample1,:);
    y = y(sample1,:);
    [n_sub,p_sub] = size(X_sub);
    
    %% Combine rankings with indices
    rankings = Z_mat(tril(true(size(Z_mat)), -1));
    [ind_i,ind_j] = find(tril(true(size(Z_mat)), -1));
    rank_tb = [rankings, ind_i, ind_j];
    
    %% Find positive/negative ranking orders
    pos_rank = rank_tb(find(rank_tb(:,1) >= 0),:);
    neg_rank = rank_tb(find(rank_tb(:,1) < 0),:);
    
    [~,pos_order] = sort(pos_rank(:,1), 'desc');
    pos_rank = pos_rank(pos_order,:);
    [~,neg_order] = sort(-neg_rank(:,1), 'desc');
    neg_rank = neg_rank(neg_order,:);
    
    %% Fowardly Find Cut-offs for Positive signals
    t_max=NaN(1,1);
    t_cut=[];
    for cut=1:size(pos_rank, 1)
        mattmp = mat_net(pos_rank(1:cut,2),pos_rank(1:cut, 3),sample1);
        Xtmp = zeros([size(mattmp, 3), 1]);
        for i = 1:size(mattmp, 3)
            Xtmp(i,:) = mean(diag(mattmp(:,:,i)));
        end
        [t_cut(cut),p_cut] = corr(Xtmp, y);
        if (abs(t_cut(cut)) > t_max) | (isnan(t_max))
            t_max = t_cut(cut);
            pvalue = p_cut;
            cut_best = cut;
        end
    end
    
    cut_best_pos = 1:cut_best;
    pos_edges = pos_rank(cut_best_pos,:);
    
     %% Fowardly Find Cut-offs for Negative signals
    t_max=NaN(1,1);
    t_cut = [];
    for cut=1:size(neg_rank, 1)
        mattmp = mat_net(neg_rank(1:cut,2),neg_rank(1:cut, 3),sample1);
        Xtmp = zeros([size(mattmp, 3), 1]);
        for i = 1:size(mattmp, 3)
            Xtmp(i,:) = mean(diag(mattmp(:,:,i)));
        end
        [t_cut(cut),p_cut] = corr(Xtmp, y);
        if (abs(t_cut(cut)) > t_max) | (isnan(t_max))
            t_max = abs(t_cut(cut));
            pvalue = p_cut;
            cut_best = cut;
        end
    end
    
    cut_best_neg = 1:cut_best;
    neg_edges = neg_rank(cut_best_neg,:);
    
end