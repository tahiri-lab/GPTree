# GPTree
## **G**enerator of **P**hylogenetic **tree**s (for supertrees)

----
> A **new** vesion of the Generator of **clusters** of phylogenetic trees with overlapping and HGT is available [here](https://github.com/tahiri-lab/GPTree/tree/TPTree_cluster).

This solution (**G**enerator of **P**hylogenetic **tree**s) generates phylogenetic trees in Newick format with a specified number of leaves and a controlled level of overlap between the trees. The generator simulates gene trees with horizontal gene transfer (HGT) and is useful for scientific experiments such as testing clustering algorithms or inferring supertrees.

If you use **GPTree** generator in your research or experiments, please consider citing the following paper:
>  Koshkarov, A., & Tahiri, N. (2023). GPTree: Generator of Phylogenetic Trees with Overlapping and Biological Events for Supertree Inference. In _BIOINFORMATICS_ (pp. 212-219).
DOI: [Link to Paper](https://doi.org/10.5220/0011697100003414)

Thank you for your contribution to the community!

The generator is based on the use of the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree).
The user has to specify several initial parameters:

*   The minimum possible number of leaves for each tree
*   The maximum possible number of leaves for each tree
*   The average level of overlap (common leaves) between the trees in the set. We will define the level of overlap between two trees as the number of common leaves (between these trees) divided by the summed length of these trees minus the number of common leaves. *Note:* for the [generator of clusters](https://github.com/tahiri-lab/GPTree/tree/TPTree_cluster) we tested another approach to calculate the average level of overlap between 2 trees.

**Initial values set by the user:**

*   **Lmin** = the minimum possible number of leaves for each tree, integer (5<=Lmin<500)
*   **Lmax** = the maximum possible number of leaves for each tree, integer (Lmin<Lmax<=500)
*   **Ngen** = the number of trees to be generated, integer (Ngen<=500)
*   **Plevel** = the average level of overlap (common leaves) between the trees in the set, in decimal notation, from 0.2 to 0.7 with steps of 0.5 (which corresponds to the range from 20% to 70%)

Currently, the generator works slow for the levels of overlap <0.2 and >0.7.

**The basic workflow:**

![The basic workflow](https://github.com/tahiri-lab/GPTree/blob/main/img/flow.png)

**Output:** the generated dataset of the specified number of trees (gene trees and/or species trees) in Newick format is saved in the folder (e.g. genetrees_50.txt file, where the number indicates the level of overlap), from which the code was launched, or in the "Files" section, if launched in colab. Examples of generated datasets see [here](https://github.com/tahiri-lab/GPTree/tree/main/test_datasets).

The Jupiter notebook also contains steps to validate the generated dataset (tree visualization, number of trees and leaves, and level of overlap).
