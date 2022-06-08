clear;
addpath '/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/code/';

% load p*p*n functional connecctome
load("/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/data/rest.mat");
% load response variable
load('/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/data/y.mat');
% load node-to-network mapping
load('/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/data/node_network_mapping.mat');

[p1,p2,n] = size(mat);
    
% Vectorize the functional connectome
X=zeros([n,p1*(p2-1)/2]);
for j=1:n
	tmp=X(:,:,j);
	X(j,:)=tmp(tril(true(size(tmp)),-1));
end
[n,p] = size(X);

% Load synthetic effects of B (large-scale)
B=load('/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/simulation/B_large.mat');
B = B.B;
B_vec = B(tril(true(size(B)),-1));
y_true=X*B_vec; 

% scNBS-lm
[scNBS_network_pos, scNBS_network_neg, scNBS_ranking, scNBS_edges_pos, scNBS_edges_neg] = scNBS(mat,y,node_network_names, node_to_network, "lm");


