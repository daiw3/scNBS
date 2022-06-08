# scNBS
The Semi-constrained Network-Based Statistic (scNBS): integrating local and global information for brain network inference

## Required Software
MATLAB: version R2020b 
MATLAB Toolbox: SparseReg; TensorReg; tensor_toolbox.

## Input data
**Functional connectome (mat)**: p*p*n, where p is the number of nodes and n is the number of samples.
**Phenotype (y)**: n*1 vector
**Node to Network mapping (node_to_network)**: a structure with 3 elements: node index (node), network index for each node (network) and network names (label).
**Edge to subnetwork mapping (node_network_names)**: a 4*[p*(p-1)/2] matrix. Each column shows the two node indices for the edge (1st and 2nd row) and corresponding two network indices (3rd and 4th row).
See more details in demo.m

## Output Results
**scNBS_pos_results**: a structure giving the marginal regression results (BETA, SD, TSTATS, P) of each subnetwork for postive associations with y. 
**scNBS_neg_results**: a structure giving the marginal regression results (BETA, SD, TSTATS, P) of each subnetwork for negative associations with y.
**scNBS_pos_edges**: a structure giving the selected edges in the functional connectome for positive associations with y. Each row consists of the two node indices for the edge.
**scNBS_neg_edges**: a structure giving the selected edges in the functional connectome for negative associations with y. Each row consists of the two node indices for the edge.
**scNBS_ranking**: a structure giving the ranking of each edge within each subnetwork. 

## Demo code for executing scNBS
See demo.m

## Simulation code
We also provided the sample simulation code we used to generate synthetic data and run scNBS (simulation folder). B_large.m is the B matrix same as the second plot of Figure 2(b) in our paper.