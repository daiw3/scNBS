clear;

map = load('data/node_network_mapping.mat');
map = map.mapping;

load('simulation/B_large.mat');

all_files = dir('results/*.mat');
nfiles = length(all_files);

scNBS_lm_pos_bonf = zeros([nfiles, 10, 10]);
scNBS_lm_neg_bonf = zeros([nfiles, 10, 10]);
scNBS_lm_pos_fdr = zeros([nfiles, 10, 10]);
scNBS_lm_neg_fdr = zeros([nfiles, 10, 10]);
scNBS_lm_pos_pval = zeros([nfiles, 10, 10]);
scNBS_lm_neg_pval = zeros([nfiles, 10, 10]);

for iter = 1:nfiles
    f1 = ['results/', all_files(iter).name];
    if(isfile(f1))
        load(f1);
    else
        continue
    end
    
    %% scNBS + Bonferroni
    
    scNBS_pos_pval = zeros(size(scNBS_network_pos));
    scNBS_neg_pval = zeros(size(scNBS_network_neg));
    for i = 1:(size(scNBS_network_pos,1))
        for j = i:size(scNBS_network_pos,2)
            if isnan(scNBS_network_pos{i,j}(4))
                scNBS_network_pos{i,j}(4)=1;
                scNBS_network_pos{i,j}(3)=0;
            end
            
            if isnan(scNBS_network_neg{i,j}(4))
                scNBS_network_neg{i,j}(4)=1;
                scNBS_network_neg{i,j}(3)=0;
            end
            
            scNBS_pos_pval(i,j) = scNBS_network_pos{i,j}(4);
            scNBS_neg_pval(i,j) = scNBS_network_neg{i,j}(4);
        end
    end
    
    scNBS_pos_bonf = zeros(size(scNBS_network_pos));
    scNBS_neg_bonf = zeros(size(scNBS_network_neg));
    tmp = pval_adjust([scNBS_pos_pval(triu(true(size(scNBS_pos_pval)))), scNBS_neg_pval(triu(true(size(scNBS_neg_pval))))], 'bonferroni') < 0.05;
    scNBS_pos_bonf(triu(true(size(scNBS_pos_pval)))) = tmp(:,1);
    scNBS_neg_bonf(triu(true(size(scNBS_neg_pval)))) = tmp(:,2);
    
    B_lm_binary = zeros(size(B));
    for i = 1:(size(scNBS_network_pos,1))
        for j = i:size(scNBS_network_pos,2)
            i_idx=find(map.network==i);
            j_idx=find(map.network==j);
            i_len=length(i_idx);
            j_len=length(j_idx);
            edges_binary = zeros([i_len, j_len]);
            if scNBS_pos_bonf(i,j) == 1
                ind = scNBS_edges_pos{i,j};
                for idx = 1:size(ind,1)
                    edges_binary(ind(idx,1), ind(idx,2))=1;
                end
            end
            if scNBS_neg_bonf(i,j) == 1
                ind = scNBS_edges_neg{i,j};
                for idx = 1:size(ind,1)
                    edges_binary(ind(idx,1), ind(idx,2))=-1;
                end
            end
            B_lm_binary(i_idx, j_idx) = edges_binary;
        end
    end
    B_lm_binary = B_lm_binary + B_lm_binary';
    B_lm_binary(B_lm_binary==-2)=-1;
    B_matrix_lm_bonf(iter,:,:) = B_lm_binary;
    scNBS_lm_pos_bonf(iter,:,:) = scNBS_pos_bonf;
    scNBS_lm_neg_bonf(iter,:,:) = scNBS_neg_bonf;
    
    
     %% scNBS + fdr
    scNBS_pos_fdr = zeros(size(scNBS_network_pos));
    scNBS_neg_fdr = zeros(size(scNBS_network_neg));
    tmp = pval_adjust([scNBS_pos_pval(triu(true(size(scNBS_pos_pval)))), scNBS_neg_pval(triu(true(size(scNBS_neg_pval))))], 'fdr') < 0.05;
    scNBS_pos_fdr(triu(true(size(scNBS_pos_pval)))) = tmp(:,1);
    scNBS_neg_fdr(triu(true(size(scNBS_neg_pval)))) = tmp(:,2);
    
    B_lm_binary = zeros(size(B));
    for i = 1:(size(scNBS_network_pos,1))
        for j = i:size(scNBS_network_pos,2)
            i_idx=find(map.network==i);
            j_idx=find(map.network==j);
            i_len=length(i_idx);
            j_len=length(j_idx);
            edges_binary = zeros([i_len, j_len]);
            if scNBS_pos_fdr(i,j) == 1
                ind = scNBS_edges_pos{i,j};
                for idx = 1:size(ind,1)
                    edges_binary(ind(idx,1), ind(idx,2))=1;
                end
            end
            if scNBS_neg_fdr(i,j) == 1
                ind = scNBS_edges_neg{i,j};
                for idx = 1:size(ind,1)
                    edges_binary(ind(idx,1), ind(idx,2))=-1;
                end
            end
            B_lm_binary(i_idx, j_idx) = edges_binary;
        end
    end
    B_lm_binary = B_lm_binary + B_lm_binary';
    B_lm_binary(B_lm_binary==-2)=-1;
    B_matrix_lm_fdr(iter,:,:) = B_lm_binary;
    scNBS_lm_pos_fdr(iter,:,:) = scNBS_pos_fdr;
    scNBS_lm_neg_fdr(iter,:,:) = scNBS_neg_fdr;
    
    scNBS_lm_pos_pval(iter,:,:) = scNBS_pos_pval;
    scNBS_lm_neg_pval(iter,:,:) = scNBS_neg_pval;

end

%% True B and B_net
Bvec= B(tril(true(size(B)),-1));
mask = tril(true(size(B_net)));

save('simu_power_large_scale_effects_results.mat', 'B_matrix_lm_bonf','B_matrix_lm_fdr', 'scNBS_lm_pos_bonf', 'scNBS_lm_neg_bonf', 'scNBS_lm_pos_fdr', 'scNBS_lm_neg_fdr', 'scNBS_lm_pos_pval', 'scNBS_lm_neg_pval', 'B', 'B_net', 'mask', 'map', 'Bvec');

