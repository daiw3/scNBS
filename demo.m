clear;
addpath '/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/code/';

% load p*p*n functional connecctome
load("/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/data/rest.mat");
% load response variable
load('/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/data/y.mat');
% load node-to-network mapping
load('/gpfs/ysm/project/scheinost/wd278/scNBS/github/scNBS/data/node_network_mapping.mat');

[scNBS_network_pos, scNBS_network_neg, scNBS_ranking, scNBS_edges_pos, scNBS_edges_neg] = scNBS(mat,y,node_network_names, node_to_network, "lm");
