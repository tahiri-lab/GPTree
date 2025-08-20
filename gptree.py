# -*- coding: utf-8 -*-
import argparse
import sys
import subprocess
import random
import pandas as pd
from ete3 import Tree
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick

def install_packages():
    """Install required packages."""
    packages = ['pandas', 'ete3', 'PyQt5', 'asymmetree']
    subprocess.check_call([sys.executable, '-m', 'pip', 'install'] + packages, stdout=subprocess.DEVNULL)

def validate_args(Lmin, Lmax, Ngen, plevel):
    """Validate user input arguments."""
    if not (5 <= Lmin < 500) or not (Lmin < Lmax <= 500) or not (3 <= Ngen <= 500):
        raise ValueError("Check if the range of the entered data is correct. Select another value, please.")
    if not (0.2 <= plevel <= 0.7):
        raise ValueError("Overlap level (plevel) must be between 0.2 and 0.7.")

def generate_initial_tree(Lmax):
    """Generate the initial species and gene tree."""
    Sf1 = te.species_tree_n(Lmax, model='yule')  # Species tree
    tree_simulator_f1 = te.GeneTreeSimulator(Sf1)
    Tf1 = tree_simulator_f1.simulate(hgt_rate=0.2, loss_rate=0.2, replace_prob=1)  # Gene tree with HGT
    ogt1 = te.prune_losses(Tf1)  # Pruned observable gene tree
    First_tree = Tree(to_newick(ogt1, reconc=False), format=1)
    First_sp1 = to_newick(Sf1, reconc=False)
    First_sp = Tree(First_sp1, format=1)  # First species tree in Newick format
    return First_tree, First_sp

def generate_tree_with_overlap(Lmin, Lmax, plevel, temp_set, tree_dataset, species_trees):
    """Generate a new gene tree with specified overlap level."""
    N1 = random.randrange(Lmin, Lmax)  # Randomly choose the number of leaves
    S1 = te.species_tree_n(N1)  # Generate next species tree
    tree_simulator = te.GeneTreeSimulator(S1)
    T1 = tree_simulator.simulate(hgt_rate=0.2, loss_rate=0.2, replace_prob=1)  # Gene tree
    ogt = te.prune_losses(T1)
    Tr1 = Tree(to_newick(ogt, reconc=False), format=1)

    overlap_tempset = []
    for i in range(len(temp_set)):
        common_leaves = set(Tr1.get_leaf_names()) & set(temp_set[i].get_leaf_names())
        overlap_level_temp = len(common_leaves) / (len(Tr1) + len(temp_set[i]) - len(common_leaves) + 0.001)
        overlap_tempset.append(overlap_level_temp)

    average_overlap = pd.Series(overlap_tempset).mean()
    if plevel - 0.01 <= average_overlap <= plevel + 0.01:
        tree_dataset.append(Tr1)
        Spt2 = Tree(to_newick(S1, reconc=False), format=1)
        species_trees.append(Spt2)
        print(f"Now we have {len(tree_dataset)} trees")

def save_dataset_to_file(filename, dataset):
    """Save the dataset of trees in Newick format to a file."""
    with open(f"{filename}.txt", "w") as file:
        for tree in dataset:
            file.write(tree.write() + '\n')

def main():
    parser = argparse.ArgumentParser(description="Generate phylogenetic trees with specified overlap.")
    parser.add_argument("Lmin", type=int, help="Minimum number of leaves (5-499)")
    parser.add_argument("Lmax", type=int, help="Maximum number of leaves (Lmin-500)")
    parser.add_argument("Ngen", type=int, help="Number of trees to generate (3-500)")
    parser.add_argument("plevel", type=float, help="Average overlap between trees (0.2-0.7)")
    args = parser.parse_args()

    try:
        validate_args(args.Lmin, args.Lmax, args.Ngen, args.plevel)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Generate the first tree
    First_tree, First_sp = generate_initial_tree(args.Lmax)

    # Initialize tree dataset
    species_trees = [First_sp]
    tree_dataset = [First_tree]
    print("Now we have 1 tree")

    # Generate trees until the desired number is reached
    while len(tree_dataset) < args.Ngen:
        generate_tree_with_overlap(args.Lmin, args.Lmax, args.plevel, tree_dataset, tree_dataset, species_trees)

    print("All done. Good job!")
    print(f"Desired number of trees = {args.Ngen}")
    print(f"Real number of trees = {len(tree_dataset)}")

    # Save generated trees
    genedataset = f"genetrees_{int(args.plevel * 100)}"
    save_dataset_to_file(genedataset, tree_dataset)

    # Save species trees if needed
    speciesdataset = f"speciestrees_{int(args.plevel * 100)}"
    save_dataset_to_file(speciesdataset, species_trees)

if __name__ == "__main__":
    install_packages()

    main()
