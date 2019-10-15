#!/usr/bin/python3

# get list of leaves of the tree
# a tree have to be rooted!
# save the tree as "all_samples_nj.tree.nex"

from Bio import Phylo
import pandas as pd
import sys

help_message = "a message"
f = input("rooted tree filename [*nex] > ")

tree_file = sys.argv[1]  # tree file name (nexus format
table_file = sys.argv[2]  # a table to be sorted

tree = Phylo.read(tree_file, "nexus")

leaves = tree.get_terminals()  # a list of Clade objects
print("Number of leaves - %i" % len(leaves))  # 528

# NB! the following reversal may be optional, see your tree
leaves = [item for item in reversed(leaves)]

# get the table
table = pd.read_csv(table_file)


