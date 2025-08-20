# -*- coding: utf-8 -*-

# Author: Aleksandr Koshkarov

import argparse
import sys
import subprocess
import random
import pandas as pd
from ete3 import Tree
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick

def install_packages():
    package = ['pandas', 'ete3', 'PyQt5', 'asymmetree']
    subprocess.check_call([sys.executable, '-m', 'pip', 'install'] + package, stdout=subprocess.DEVNULL)

def validate_args(k, Lmin, Lmax, Ngen, plevel):
    if not (1 <= k <= 100) or not (5 <= Lmin < Lmax <= 500) or not (3 <= Ngen <= 500):
        raise ValueError("Invalid parameter values. Please check the range and constraints.")
    if not (0.2 <= plevel <= 0.7):
        raise ValueError("Overlap level (plevel) must be between 0.2 and 0.7.")

# Function to rename leaves of a tree with "L" prefix
def rename_leaves(tree):
    for leaf in tree:
        leaf.name = f"L{leaf.name}"  # Add "L" before each leaf name
    return tree

def gptree_speciestree(Lmin, Lmax):
    N1 = random.randrange(Lmin, Lmax)  # Randomly choose the number of leaves for the generated tree
    return te.species_tree_n(N1)

# Updated function to generate one gene tree and rename leaves
def gptree_genetree(S1, hgt_rate=0.2, loss_rate=0.2, replace_prob=0.9):
    tree_simulator = te.GeneTreeSimulator(S1)
    T1 = tree_simulator.simulate(hgt_rate=hgt_rate, loss_rate=loss_rate, replace_prob=replace_prob)
    ogt = te.prune_losses(T1)
    Gn_Tree = Tree(to_newick(ogt, reconc=False), format=1)
    return rename_leaves(Gn_Tree)

def calculate_overlap(tree1, tree2):
    common_leaves = set(tree1.get_leaf_names()) & set(tree2.get_leaf_names())
    overlap_level = ((len(common_leaves) / len(tree1)) + (len(common_leaves) / len(tree2))) / 2 # Sørensen–Dice
    #overlap_level = len(common) / (len(tree1) + len(tree2) - len(common)) # uncomment for the Jaccard coefficient
     return overlap_level

def gptree_cluster_gene(sptree, Ngen, plevel):
    cluster_dataset = [gptree_genetree(sptree)]  # Add the first gene tree
    print("Now we have 1 tree")
    
    while len(cluster_dataset) < Ngen:
        gene_tree_next = gptree_genetree(sptree)
        overlap_levels = [calculate_overlap(gene_tree_next, gt) for gt in cluster_dataset]
        average_overlap = sum(overlap_levels) / len(overlap_levels)
        
        if plevel - 0.01 <= average_overlap <= plevel + 0.01:
            cluster_dataset.append(gene_tree_next)
            print(f"Now we have {len(cluster_dataset)} trees")
    
    return cluster_dataset

def main():
    parser = argparse.ArgumentParser(description="Generate clusters of phylogenetic trees with specified overlap.")
    parser.add_argument("k", type=int, help="Number of clusters (1-100)")
    parser.add_argument("Lmin", type=int, help="Minimum number of leaves (5-499)")
    parser.add_argument("Lmax", type=int, help="Maximum number of leaves (Lmin-500)")
    parser.add_argument("Ngen", type=int, help="Number of trees per cluster (3-500)")
    parser.add_argument("plevel", type=float, help="Average overlap between trees (0.2-0.7)")
    
    args = parser.parse_args()
    
    try:
        validate_args(args.k, args.Lmin, args.Lmax, args.Ngen, args.plevel)
    except ValueError as e:
        print(e)
        sys.exit(1)
    
    tree_cluster_dataset = f"trees_{args.k}_{args.Lmin}_{args.Lmax}_{args.Ngen}_{int(args.plevel * 100)}"
    with open(f"{tree_cluster_dataset}.txt", "w") as file:
        for i in range(1, args.k + 1):
            print(f"For cluster {i}")
            species_tree_k = gptree_speciestree(args.Lmin, args.Lmax)
            cluster_k = gptree_cluster_gene(species_tree_k, args.Ngen, args.plevel)
            
            for tree in cluster_k:
                file.write(tree.write() + '\n')
    
    print(f"All done. Generated trees saved to {tree_cluster_dataset}.txt")

if __name__ == "__main__":
    install_packages()
    main()

