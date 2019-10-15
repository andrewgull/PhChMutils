#!/usr/bin/env python3

# get list of leaves of the tree
# a tree have to be rooted!
# save the tree as "all_samples_nj.tree.nex"

from Bio import Phylo
import pandas as pd
import sys

help_message = "A script for reordering table columns according to a tree\n" \
               "Arguments:\n" \
               "1 - rooted tree in nexus format\n" \
               "2 - SNP table, comma or TAB separated file with '.csv' extension\n"

if len(sys.argv) < 3:
    print(help_message)
    sys.exit()

tree_file = sys.argv[1]  # tree file name (nexus format
table_file = sys.argv[2]  # a table to be sorted

tree = Phylo.read(tree_file, "nexus")

leaves = tree.get_terminals()  # a list of Clade objects
print("Number of leaves - %i" % len(leaves))

# NB! the following reversal may be optional, see your tree
leaves = [item for item in reversed(leaves)]

print("%s ... %s" % (leaves[1], leaves[-1]))

# get the table
table = pd.read_csv(table_file, sep=None, engine='python')

table = table.reindex(list(table.columns)[1:13] + leaves, axis=1)

table.to_csv(table_file[:-4] + "_treeorder.csv", sep="\t", na_rep="", index=False)

print("Done!")
