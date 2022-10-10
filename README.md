# GPTreeCluster
## Generator of clusters of phylogenetic trees with overlapping and HGT

This generator can be used to generate a specified number of clusters of phylogenetic trees in Newick format with a *variable number of leaves*  and with some level of *overlap* between trees in clusters. With this tool, the user can generate a dataset with clusters of gene trees (particularly, gene trees with **horizontal gene transfer** implemented), which is saved in txt, with the possibility of its further use in their scientific experiments (e.g., testing classification algorithms or inference supertrees).

The generator is based on the use of the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree).
The user has to specify several initial parameters:

*   The number of clusters.
*   The minimum possible number of leaves for trees in a cluster.
*   The maximum possible number of leaves for trees in a cluster.
*   The average level of overlap (common leaves) between the trees in each cluster. We will define the level of overlap between trees $T_1$ and $T_2$ as the average between 2 values: (1) the number of common leaves in trees $T_1$ and $T_2$, divided by the number of leaves in tree $T_1$, (2)  the number of common leaves in trees $T_1$ and $T_2$, divided by the number of leaves in tree $T_2$.

**Important note:** For calculation the level of overlapping we can use several metrics such as the Jaccard Similarity, also called the Jaccard Index or Jaccard Similarity Coefficient (given two sets, A and B, the Jaccard Similarity is defined as the size of the intersection of set A and set B (i.e. the number of common elements) over the size of the union of set A and set B (i.e. the number of unique elements)) or the Overlap Coefficient, also known as the Szymkiewicz–Simpson coefficient, is defined as the size of the union of set A and set B over the size of the smaller set between A and B. To get more stable results we calculate (in this version) the mean overlap using both sets.

**Initial values set by the user:**

*   **k** = the number of clusters, integer (1<=k<100).
*   **Lmin** = the minimum possible number of leaves for trees in a cluster, integer (5<=Lmin<500)
*   **Lmax** = the maximum possible number of leaves for trees in a cluster, integer (Lmin<Lmax<=500)
*   **Ngen** = the number of trees to be generated for each cluster, integer (Ngen<=500)
*   **Plevel** = the average level of overlap (common leaves) between trees in a cluster, in decimal notation, from 0.2 to 0.7 with steps of 0.5 (which corresponds to the range from 20% to 70%)

Currently, the generator works very slow for the levels of overlap <0.2 and >0.7.

**The basic workflow:**

![The basic workflow](https://github.com/tahiri-lab/GPTree/blob/TPTree_cluster/img/flow.png)

**Output:** The generated dataset of the specified number of trees (separated by clusters) in Newick format is saved in the folder (e.g. trees_3_6_50.txt file, where the numbers indicate the number of clusters, number of trees in each cluster, and level of overlapping).
