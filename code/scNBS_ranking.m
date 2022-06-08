function [Z_score] = scNBS_ranking(mat,y,sample1, colnames, node_to_network, r_method)
    
    n = size(mat, 3);
    n_networks=histc(node_to_network.network,unique(node_to_network.network));
        
    Z_score = cell(length(n_networks), length(n_networks));
    for net1 = 1:length(n_networks)
        for net2 = net1:length(n_networks)
            node_idx1= node_to_network.node(find(node_to_network.network == net1));
            node_idx2= node_to_network.node(find(node_to_network.network == net2));
            X_subnetwork = mat(node_idx1, node_idx2,:);
            
            if length(node_idx1) == length(node_idx2)
                X_subnetwork_vec = zeros([n,length(node_idx1)*(length(node_idx1)-1)/2]);
                for j=1:n
                    tmp=X_subnetwork(:,:,j);
                    X_subnetwork_vec(j,:)=tmp(tril(true(size(tmp)),-1));
                end 
            else
                X_subnetwork_vec = zeros([n,length(node_idx1)*length(node_idx2)]);
                for j=1:n
                    tmp=X_subnetwork(:,:,j);
                    X_subnetwork_vec(j,:)=tmp(:);
                end 
            end
            
            %% perform ranking, depending on different methods
            
            % OLS
            if r_method=="lm"
                [Zvalue,pvalue] = corr(X_subnetwork_vec(sample1,:),y(sample1,:));
                %Z_score{net1,net2} = Zvalue;
                clear pvalue;
                %clear Zvalue;
                %disp('Done.\n')
            
            % Ridge
            elseif r_method=="ridge"
                [B,FitInfo] = lasso(X_subnetwork_vec(sample1,:),y(sample1,:),'Alpha',1e-9,'CV',10);
                idxLambda1SE = FitInfo.Index1SE;
                Zvalue = B(:,idxLambda1SE);
            
            % Lasso
            elseif r_method=="lasso"
                [B,FitInfo] = lasso(X_subnetwork_vec(sample1,:),y(sample1,:),'Alpha',1,'CV',10);
                idxLambda1SE = FitInfo.Index1SE;
                Zvalue = B(:,idxLambda1SE);
            
            % Matrix 
            elseif r_method=="matrix"
                B = cell(3,1);
                BIC = zeros(3,1);
                %tic;
                for r = 1:3
                  [~,beta_rk_tmp,glmstats_tmp,~] = kruskal_reg(ones(length(sample1),1),X_subnetwork(:,:,sample1),y(sample1,:),r,'normal');
                  BIC(r) = glmstats_tmp{end}.BIC;
                  B{r} = beta_rk_tmp;
                end
                %toc;
                [~,bic_ind] = min(BIC);
                Zvalue = double(B{bic_ind});
                clear beta_rk_tmp;
                clear glmstats_tmp;
                clear B;
                
                %% covert 2D z_score to 1D
                if length(node_idx1) == length(node_idx2)
                    tmp = (Zvalue+Zvalue')/2;
                    Ztmp = tmp(tril(true(size(tmp)),-1));
                else
                    Ztmp = Zvalue(:);
                end
                Zvalue = Ztmp;
                
            
            % Matrix Lasso
            elseif r_method=="matrix_lasso"
                B = cell(3,1);
                BIC = zeros(3,1);
                ranks = 1:3;
                for r = ranks
                  [~,beta_rk_tmp,glmstats_tmp,~] = kruskal_reg(ones(length(sample1),1),X_subnetwork(:,:,sample1),y(sample1,:),r,'normal');
                  BIC(r) = glmstats_tmp{end}.BIC;
                  B{r} = beta_rk_tmp;
                end
                [~,bic_ind] = min(BIC);
                Btmp = B{bic_ind};
                clear beta_rk_tmp;
                clear glmstats_tmp;
                clear B;
                
                % Set lasso penalty and tuning parameter values
                pentype = 'enet';
                penparam = 1;
                [~,~,stats] = matrix_sparsereg(ones(length(sample1),1),X_subnetwork(:,:,sample1),y(sample1,:),inf,'normal');
                maxlambda = stats.maxlambda*.95;
                %maxlambda = 1000;
                gridpts = 10;
                lambdas = zeros(1,gridpts);
                gs = 2/(1+sqrt(5));
                B = cell(1,gridpts);
                BIC = zeros(1,gridpts);
                %tic;
                for i=1:gridpts
                    lambda = maxlambda*gs^(i-1);
                    lambdas(i) = lambda;
                    [~,B{i},~, stats] = kruskal_sparsereg(ones(length(sample1),1),X_subnetwork(:,:,sample1),y(sample1,:),ranks(bic_ind),'normal',...
                                               lambda,pentype,penparam,'B0', Btmp);
                    BIC(i) = stats{end}.BIC;
                end
                %toc;
                [~,bic_ind] = min(BIC);
                Zvalue = double(B{bic_ind});
                
                %% covert 2D z_score to 1D
                if length(node_idx1) == length(node_idx2)
                    tmp = (Zvalue+Zvalue')/2;
                    Ztmp = tmp(tril(true(size(tmp)),-1));
                else
                    Ztmp = Zvalue(:);
                end
                Zvalue = Ztmp;
            end
            
            
            
            
            %% covert 1D z_score to 2D
            if length(node_idx1) == length(node_idx2)
               tmp=zeros([length(node_idx1), length(node_idx1)]);
               ind = find(tril(ones(length(node_idx1),length(node_idx1)),-1));
               tmp(ind) = Zvalue;
               Z_mat = tmp + tmp';
            else
               Z_mat = reshape(Zvalue, [length(node_idx1), length(node_idx2)]);
            end
            Z_score{net1,net2} = Z_mat;
            % Z_score{net1,net2} = Zvalue;
            
        end
    end
end