## Generator of clusters of phylogenetic trees with different number of overlapping leaves

To run:

`python gptreecluster.py k Lmin Lmax Ngen p `

where:

*k* - number of clusters  \
*Lmin* - minimum possible number of leaves for trees in a cluster \
*Lmax* - maximum possible number of leaves for trees in a cluster \
*Ngen* - number of tree to generate (in each cluster) \
*p* - average level of overlap (common leaves) between the trees in each cluster

The recommended intervals of possible values are as follows:

 * $1 \leq k \leq 100$, integer; 
 * $5 \leq L_{min}<500$, integer; 
 * $L_{min} < L_{max} \leq 500$, integer; 
 * $N_{gen} \leq 500$, integer; 
 * $0.2 \leq p \leq 0.7$, float.

Example:

`python gptreecluster.py 3 15 25 40 0.5 `

The generated dataset of the specified number of trees (separated by clusters) in Newick format is saved in the folder (e.g. trees_3_6_50.txt file, where the numbers indicate the number of clusters, number of trees in each cluster, and level of overlapping), from which the code was launched.
