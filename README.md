# GPTree
## **G**enerator of **P**hylogenetic **tree**s (for supertrees)

----
> :bulb: A **new** vesion of the Generator of **clusters** of phylogenetic trees with overlapping and HGT is available [here](https://github.com/tahiri-lab/GPTree/tree/TPTree_cluster).

This solution (**G**enerator of **P**hylogenetic **tree**s) generates phylogenetic trees in Newick format with a specified number of leaves and a controlled level of overlap between the trees. The generator simulates gene trees with horizontal gene transfer (HGT) and is useful for scientific experiments such as testing clustering algorithms or inferring supertrees.

:bookmark_tabs: If you use **GPTree** generator in your research or experiments, please consider citing the following paper:
>  Koshkarov, A., & Tahiri, N. (2023). GPTree: Generator of Phylogenetic Trees with Overlapping and Biological Events for Supertree Inference. In _BIOINFORMATICS_ (pp. 212-219).
DOI: [Link to Paper](https://doi.org/10.5220/0011697100003414)

:trophy: Thank you for your contribution to the community!

The generator is based on the use of the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree).

## Features

- Generates phylogenetic trees with **horizontal gene transfer (HGT)**.
- Allows users to specify the number of leaves and overlap level between trees.
- Outputs gene trees and species trees in **Newick format**.
- Designed to handle large datasets with configurable parameters.

## Requirements

The script depends on the following Python libraries:
- `ete3`
- `PyQt5`
- `asymmetree`
- `pandas`

## Input Parameters

The user needs to provide several initial parameters:

1. **Lmin**: Minimum number of leaves per tree (integer, 5 ≤ Lmin < 500).
2. **Lmax**: Maximum number of leaves per tree (integer, Lmin < Lmax ≤ 500).
3. **Ngen**: Number of trees to generate (integer, 3 ≤ Ngen ≤ 500).
4. **Plevel**: Average level of overlap (common leaves) between trees, as a decimal (0.2 ≤ plevel ≤ 0.7).

The overlap level between trees is calculated based on the number of common leaves between them, with additional controls to ensure the desired level of overlap.

Currently, the generator works slow for the levels of overlap <0.2 and >0.7.

**The basic workflow:**

![The basic workflow](https://github.com/tahiri-lab/GPTree/blob/main/img/flow.png)

Examples of generated datasets see [here](https://github.com/tahiri-lab/GPTree/tree/main/test_datasets).

The [Jupiter notebook](https://github.com/tahiri-lab/GPTree/blob/main/Overlap_Phyltree_generator_for_supertree_inference.ipynb) also contains steps to validate the generated dataset (tree visualization, number of trees and leaves, and level of overlap).
