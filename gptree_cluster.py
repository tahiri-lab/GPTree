# -*- coding: utf-8 -*-

# Author: Aleksandr Koshkarov, UdeS

import argparse
import sys
import subprocess
import random
import pandas as pd  # (optional) unused
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
    for leaf in tree.iter_leaves():  # iterate over leaves only
        leaf.name = f"L{leaf.name}"
    return tree

def gptree_speciestree(Lmin, Lmax):
    N1 = random.randrange(Lmin, Lmax)
    return te.species_tree_n(N1)

# Function to generate one gene tree and rename leaves
def gptree_genetree(S1, hgt_rate=0.2, loss_rate=0.2, replace_prob=0.9):
    tree_simulator = te.GeneTreeSimulator(S1)
    T1 = tree_simulator.simulate(hgt_rate=hgt_rate, loss_rate=loss_rate, replace_prob=replace_prob)
    ogt = te.prune_losses(T1)
    Gn_Tree = Tree(to_newick(ogt, reconc=False), format=1)
    return rename_leaves(Gn_Tree)

def calculate_overlap(tree1, tree2):
    # Jaccard over leaf sets
    leaves1 = set(tree1.get_leaf_names())
    leaves2 = set(tree2.get_leaf_names())
    common = leaves1 & leaves2
    union = leaves1 | leaves2
    return len(common) / len(union) if union else 0.0

def gptree_cluster_gene(sptree, Ngen, plevel):
    cluster_dataset = [gptree_genetree(sptree)]
    print("Now we have 1 tree")
    while len(cluster_dataset) < Ngen:
        gene_tree_next = gptree_genetree(sptree)
        overlap_levels = [calculate_overlap(gene_tree_next, gt) for gt in cluster_dataset]
        average_overlap = sum(overlap_levels) / len(overlap_levels)
        if plevel - 0.01 <= average_overlap <= plevel + 0.01:
            cluster_dataset.append(gene_tree_next)
            print(f"Now we have {len(cluster_dataset)} trees")
    return cluster_dataset

# Topology signature helpers

def tree_topology_signature(ete_tree: Tree):
    """
    Return a hashable topology signature:
    - leafset: frozenset of leaf names
    - splits: frozenset of frozensets, each being the set of leaf names under an internal node
    Branch lengths are ignored.
    """
    leafset = frozenset(ete_tree.get_leaf_names())
    splits = set()
    for node in ete_tree.traverse("postorder"):
        if not node.is_leaf():
            clade = frozenset(leaf.name for leaf in node.iter_leaves())
            if 1 < len(clade) < len(leafset):
                splits.add(clade)
    return (leafset, frozenset(splits))

def species_topology_signature(species_tree_obj):
    """
    Convert the asymmetree species tree to Newick, parse with ETE, and build the signature.
    """
    nwk = to_newick(species_tree_obj, reconc=False)
    ete_t = Tree(nwk, format=1)
    return tree_topology_signature(ete_t)

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

    # Keep track of seen species-tree topologies
    seen_species_topologies = set()

    tree_cluster_dataset = f"trees_{args.k}_{args.Lmin}_{args.Lmax}_{args.Ngen}_{int(args.plevel * 100)}"
    with open(f"{tree_cluster_dataset}.txt", "w") as file:
        for i in range(1, args.k + 1):
            print(f"For cluster {i}")

            # Keep generating until we find a species tree with a new topology
            while True:
                species_tree_k = gptree_speciestree(args.Lmin, args.Lmax)
                sig = species_topology_signature(species_tree_k)
                if sig in seen_species_topologies:
                    continue  # duplicate; try again
                seen_species_topologies.add(sig)
                break

            cluster_k = gptree_cluster_gene(species_tree_k, args.Ngen, args.plevel)

            for tree in cluster_k:
                file.write(tree.write() + '\n')

    print(f"All done. Generated trees saved to {tree_cluster_dataset}.txt")

if __name__ == "__main__":
    install_packages()
    main()
