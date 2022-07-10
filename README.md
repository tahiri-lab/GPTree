# GPTree
## Phylogenetic tree generator (for supertrees)

## Générateur d'arbres phylogénétiques (pour les super-arbres)

This [generator](https://github.com/tahiri-lab/GPTree/blob/main/Overlap_Phyltree_generator_for_supertree_inference.ipynb) can be used to generate a specified number of phylogenetic trees in Newick format with a variable number of leaves and with some level of overlap between trees. With this tool, the user can generate a dataset with trees (particularly, gene trees with **[horizontal gene transfer](https://github.com/tahiri-lab/GPTree/tree/main/HGT_test)** implemented), which is saved in txt, with the possibility of its further use in their scientific experiments (e.g., testing classification algorithms or inference supertrees).

The generator is based on the use of the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree).
The user has to specify several initial parameters:

*   The minimum possible number of leaves for each tree
*   The maximum possible number of leaves for each tree
*   The average level of overlap (common leaves) between the trees in the set. We will define the level of overlap between two trees as the number of common leaves (between these trees) divided by the summed length of these trees minus the number of common leaves.

**Initial values set by the user:**

*   **Lmin** = the minimum possible number of leaves for each tree, integer (5<=Lmin<500)
*   **Lmax** = the maximum possible number of leaves for each tree, integer (Lmin<Lmax<=500)
*   **Ngen** = the number of trees to be generated, integer (Ngen<=500)
*   **Plevel** = the average level of overlap (common leaves) between the trees in the set, in decimal notation, from 0.2 to 0.7 with steps of 0.5 (which corresponds to the range from 20% to 70%)

Currently, the generator works very slow for the levels of overlap <0.2 and >0.7.

**The basic workflow:**

![The basic workflow](https://github.com/tahiri-lab/GPTree/blob/main/img/flow.png)

**Output:** the generated dataset of the specified number of trees (gene trees and/or species trees) in Newick format is saved in the folder (e.g. genetrees_50.txt file, where the number indicates the level of overlap), from which the code was launched, or in the "Files" section, if launched in colab. Examples of generated datasets see [here](https://github.com/tahiri-lab/GPTree/tree/main/test_datasets).

The Jupiter notebook also contains steps for checking the generated dataset (tree visualization, number of trees and leaves, level of overlap).
