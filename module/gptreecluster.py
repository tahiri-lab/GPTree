# -*- coding: utf-8 -*-
"""
##Generator of clusters of phylogenetic trees with overlapping and HGT

This generator can be used to generate a specified number of clusters of phylogenetic trees in Newick format with a *variable number of leaves*  and with some level of *overlap* between trees in clusters. With this tool, the user can generate a dataset with clusters of gene trees (particularly, gene trees with **horizontal gene transfer** implemented), which is saved in txt, with the possibility of its further use in their scientific experiments (e.g., testing classification algorithms or inference supertrees).

Contributors: AK and NT

The generator is based on the use of the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree).
The user has to specify several initial parameters:

*   The number of clusters.
*   The minimum possible number of leaves for trees in a cluster.
*   The maximum possible number of leaves for trees in a cluster.
*   The average level of overlap (common leaves) between the trees in each cluster. We will define the level of overlap between trees $T_1$ and $T_2$ as the average between 2 values: (1) the number of common leaves in trees $T_1$ and $T_2$, divided by the number of leaves in tree $T_1$, (2)  the number of common leaves in trees $T_1$ and $T_2$, divided by the number of leaves in tree $T_2$.

To generate species and gene trees with horizintal gene transfer we use this library: https://github.com/david-schaller/AsymmeTree

**Generator**

"""
import argparse
import sys
import subprocess

# implement pip as a subprocess:
package = ['pandas', 'ete3', 'PyQt5', 'asymmetree']
subprocess.check_call([sys.executable, '-m', 'pip', 'install'] + package, stdout=subprocess.DEVNULL)

# process output with an API in the subprocess module:
reqs = subprocess.check_output([sys.executable, '-m', 'pip'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
#print(installed_packages)


#Importing the required libraries
import PyQt5
from ete3 import Tree, PhyloTree, TreeStyle
import random
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick
import pandas as pd

"""**Setting the required initial values by the user**

*   **k** = the number of clusters, integer (1<=k<100).

*   **Lmin** = the minimum possible number of leaves for trees in a cluster, integer (5<=Lmin<500)
*   **Lmax** = the maximum possible number of leaves for trees in a cluster, integer (Lmin<Lmax<=500)
*   **Ngen** = the number of trees to be generated for each cluster, integer (Ngen<=500)
*   **Plevel** = the average level of overlap (common leaves) between trees in a cluster, in decimal notation, from 0.2 to 0.7 with steps of 0.5 (which corresponds to the range from 20% to 70%)

Currently, the generator works very slow for the levels of overlap <0.2 and >0.7.
"""

parser=argparse.ArgumentParser(description="Use these arguments: gptreecluster.py k Lmin Lmax Ngen p")
#args = parser.parse_args()

# sys.argv includes a list of elements starting with the program
if len(sys.argv) < 5:
    parser.print_usage()
    print("""
    Instruction:
    To generate clusters, the user has to specify the following parameters:
      gptreecluster.py k Lmin Lmax Ngen p
    -- k - number of clusters;
    -- Lmin - minimum possible number of leaves for trees in a cluster;
    -- Lmax - maximum possible number of leaves for trees in a cluster;
    -- Ngen - number of tree to generate (in each cluster);
    -- p - average level of overlap (common leaves) between the trees in each cluster.
    The generated dataset of the specified number of trees (separated by clusters) 
    in Newick format is saved in the folder (e.g. trees_3_40_50.txt file, where 
    the numbers indicate the number of clusters, number of trees in each cluster, 
    and level of overlapping), from which the code was launched.
    """)
    sys.exit(1)

k = int(sys.argv[1])
Lmin = int(sys.argv[2])
Lmax = int(sys.argv[3])
Ngen = int(sys.argv[4])
plevel = float(sys.argv[5])


#Checking values
if k < 1 or k > 100 or Lmin < 5 or Lmin >499 or Lmin > Lmax or Lmin==Lmax or Lmax >500 or Ngen < 3 or Ngen >500:
  print('Check if the range of the entered data is correct. Select another value, please')
  raise ValueError

# a function to generate a species tree with random number of leaves

def gptree_speciestree():
  N1 = random.randrange(Lmin, Lmax) #we randomly chose the number of leaves for the generated tree in the specified range
  S1 = te.species_tree_n(N1)
  return S1

# a function to generate one gene tree based on a species tree

def gptree_genetree(S1):
  tree_simulator = te.GeneTreeSimulator(S1)
  T1 = tree_simulator.simulate(hgt_rate=0.2, loss_rate=0.2,  replace_prob=0.9) # we set the horizontal gene transfer rate to 0.2
  ogt = te.prune_losses(T1)
  ogt = to_newick(ogt, reconc=False)
  Gn_Tree = Tree(ogt, format=1)
  return Gn_Tree

# a function to generate a set of gene trees based on a species tree

def gptree_cluster_gene(sptree):
  
  cluster_dataset = [] #declare an empty list for gene trees
  gene_tree1 = gptree_genetree(sptree)
  cluster_dataset.append(gene_tree1) #add the first gene tree
  common_leaves_temp = []
  overlap_set1 = []
  while len(cluster_dataset) < Ngen:
    temp_set = cluster_dataset #temp_set is to store trees temporarily
    overlap_tempset = overlap_set1
    gene_tree_next = gptree_genetree(sptree)
    # in the next step we will check the level of overlap and add trees to the intermediate set
    overlap_tempset =[]
    common_leaves_temp = []
    for i in range (0, len(temp_set)):
      rf, max_rf, common_leaves, parts_t1, parts_t2, discard_t1, discart_t2 = gene_tree_next.robinson_foulds(temp_set[i], unrooted_trees=True)
      # we used the robinson_foulds function to get common leaves for each pairs of trees
      overlap_level_temp = len(common_leaves)/(len(gene_tree_next)+len(temp_set[i]) - len(common_leaves)) # overlap level as the Jaccard Similarity Coefficient
      #overlap_level_temp = len(common_leaves)/(min(len(gene_tree_next), len(temp_set[i]))) #overlap level as the Overlap Coefficient
      #overlap_level_temp = ((len(common_leaves)/len(gene_tree_next))+(len(common_leaves)/len(temp_set[i])))/2 # overlap level considering average
      overlap_tempset.append(overlap_level_temp)
    average_overlap = pd.Series(overlap_tempset) #calculation of the average level of overlap in the current version of our dataset
    if (average_overlap.mean() <= plevel+0.01) and (average_overlap.mean() >= plevel-0.01):
      overlap_set1.append(overlap_level_temp)
      cluster_dataset.append(gene_tree_next) #we add the generated gene tree to our dataset if this tree satisfies our conditions
      print("Now we have ", len(cluster_dataset), " trees") #we display the current number of trees in our dataset to see the progress
  #print("All done. Good job!")
  print("Now we have ", len(cluster_dataset), " trees")
  return cluster_dataset

# Generate k clusters of trees and save them to the txt file

tree_cluster_dataset = "trees_%s_%s_%s_%s_%s" % (k, Lmin, Lmax, Ngen, int(plevel*100))
with open(r"%s.txt" % tree_cluster_dataset, "w") as file:
  for i in range(1, k+1):
    print("For cluster ", i)
    species_tree_k = gptree_speciestree()
    cluster_k = gptree_cluster_gene(species_tree_k)
  # now we write the generated gene trees to a txt file (trees in Newick format).
    #file.write("Cluster # %s" % (i)  + '\n') #if you don't need to separate clusters by the word "Cluster #", you can comment this line

  #write trees into the file line by line
    for tree in cluster_k:
      file.write(tree.write() + '\n')
file.close() #close the file for writing
print("All done. Good job!")
"""The generated dataset of the specified number of trees (separated by clusters) in Newick format is saved in the folder (e.g. trees_3_6_50.txt file, where the numbers indicate the number of clusters, number of trees in each cluster, and level of overlapping), from which the code was launched.

"""
