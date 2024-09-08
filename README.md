# GPTreeCluster

## Generator of Clusters of Phylogenetic Trees with Overlapping and Horizontal Gene Transfer (HGT)

This tool generates clusters of phylogenetic trees in Newick format, with a *variable number of leaves* and a *configurable overlap* between the trees in each cluster. The trees also incorporate **horizontal gene transfer (HGT)**. The generated dataset can be used for scientific experiments such as testing classification algorithms or inferring supertrees.

### Features:
- Generates phylogenetic trees with **HGT**.
- Allows users to specify the number of clusters and trees, as well as the overlap level between trees.
- Saves output in Newick format for easy integration with other tools.
- Based on the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree).

### Prerequisites

The following Python libraries are required:
- `pandas`
- `ete3`
- `PyQt5`
- `asymmetree`

These can be installed using pip:
```
pip install pandas ete3 PyQt5 asymmetree
```

### Input Parameters:

1. **k**: Number of clusters (integer, 1 ≤ k ≤ 100).
2. **Lmin**: Minimum number of leaves per tree (integer, 5 ≤ Lmin < 500).
3. **Lmax**: Maximum number of leaves per tree (integer, Lmin < Lmax ≤ 500).
4. **Ngen**: Number of trees to generate per cluster (integer, Ngen ≤ 500).
5. **p**: Average level of overlap (common leaves) between trees in a cluster (float, 0.2 ≤ p ≤ 0.7).

### Overlap Calculation:

The level of overlap between trees can be calculated using different metrics, such as:

- **Jaccard Similarity (default)**: Measures the ratio of common leaves between two trees to the total number of unique leaves.
  
  $$Jaccard(A, B) = \frac{|A \cap B|}{|A \cup B|}$$

- **Overlap Coefficient**: Measures the ratio of common leaves relative to the smaller set of leaves between two trees.
  
  $$Overlap(A, B) = \frac{|A \cap B|}{\min(|A|, |B|)}$$

- **Mean Overlap**: The average of the common leaves relative to both trees.
  
  $$MeanOverlap(A, B) = \frac{\frac{|A \cap B|}{|A|} + \frac{|A \cap B|}{|B|}}{2}$$

### Usage:

To run the generator, use the following command:

```bash
python gptree_cluster.py k Lmin Lmax Ngen p
```

Where:
- `k` is the number of clusters,
- `Lmin` is the minimum number of leaves per tree,
- `Lmax` is the maximum number of leaves per tree,
- `Ngen` is the number of trees per cluster,
- `p` is the average level of overlap between the trees in each cluster.

### Example:

To generate 3 clusters with trees having between 15 and 25 leaves, 40 trees per cluster, and an overlap level of 50%:

```bash
python gptree_cluster.py 3 15 25 40 0.5
```

This will output a file named `trees_3_15_25_40_50.txt`, where the numbers represent the input parameters.

### Output:

The generated dataset is saved in a text file in Newick format, where each line represents a tree. The filename follows the format: 

```
trees_<k>_<Lmin>_<Lmax>_<Ngen>_<p>.txt
```

For example, a file named `trees_3_15_25_40_50.txt` would indicate:
- **3 clusters**
- **15–25 leaves per tree**
- **40 trees per cluster**
- **50% average overlap between trees**

### Workflow Overview:

Here is a high-level workflow of the generator:

![Workflow](https://github.com/tahiri-lab/GPTree/blob/GPTreeCluster/img/flow.png)

### Notes:

- The generator may slow down for overlap levels below 0.2 or above 0.7.
- Ensure that the input parameters fall within the specified ranges for optimal performance.
