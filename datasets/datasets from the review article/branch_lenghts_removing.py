# Branch lenghts removing

# In order to use RF(-) and RF(+) we need phylogenetic trees without branch lengths and in Newick format. When we use GPTree Cluster generator, the generated data includes branch lengths.

# Here we remove branch lengths and bootstrap values, and rename leaf names by adding "L" before each leaf name

import re

# Read in the input file
with open('clusters3.txt', 'r') as input_file:
    input_trees = input_file.readlines()

# Remove branch lengths and bootstrap values from each tree, and add "L" before each leaf name
output_trees = []
for tree in input_trees:
    # Use regular expressions to remove branch lengths and bootstrap values, and add "L" before each leaf name
    output_tree = re.sub(r':\d+(\.\d+)?', '', tree)  # remove branch lengths
    output_tree = re.sub(r'\)[\d\.]+', ')', output_tree)  # remove bootstrap values
    output_tree = re.sub(r'([^(,:]+)', r"L\1", output_tree)  # add "L" before each leaf name
    output_trees.append(output_tree)

# Write the modified trees to a new file
with open('trees_formatted.txt', 'w') as output_file:
    output_file.writelines(output_trees)