#!/usr/bin/python3

# The script is for assigning each sample to phylo lineage by specific SNPs
# Also it removes samples that got different lineages' signatures
# and outputs it to txt file
# possible args: phylo_table, ignore list, number of characters to cut off from sample name (here 11)

import pandas as pd
import glob
from collections import defaultdict
import re
import numpy as np
import sys


def dict_reverse(d):
    """
    function for reversing the dict of lists
    :param d: dict to reverse
    :return: reversed dict
    """
    d_rev = defaultdict(list)
    for k in d.keys():
        for v in d[k]:
            d_rev[v].append(k)
    return d_rev

def filter_mixed_samples(smpl_dict, *params):
    """
    function for removing mixed samples from samples dict
    adding them to special list
    :param smpl_dict: sample dict - {sample1:[lineage1, lineage2]}
    :param params: lineages names to be ignored
    :return: cleaned sample dict
    """
    bad_samples = []
    ignore_list = params
    # ignore = ['lineage4', 'lineage4.9']
    for key in smpl_dict.keys():
        # first remove the ignores
        lst = [item for item in smpl_dict[key] if item not in ignore_list]
        #  if each key has more than one value
        if len(lst) > 1:
            # leave only 'lineageX' names of lineages
            lst = [re.search(r'(lineage[0-9])', item) for item in lst]
            lst = [item.group(0) for item in lst]
        # check if tems are NOT identical
        # if true put them to separate list
        if lst[1:] != lst[:-1]:
            bad_samples.append(key)

    # after iteration remove bad samples from the dict
    for bad in bad_samples:
        # use del not to generate None values
        del smpl_dict[bad]
        # sample_dict.pop(bad, None)
    return smpl_dict, bad_samples

def make_phylo_table(ph_table, filter):
    # read the lineage table
    LS_snps = pd.read_csv(ph_table)

    # collect VCF file names
    files = glob.glob("*.vcf")

    # get pos, lineages' names and allele changes
    snp_positions = LS_snps["Position"]
    lins = LS_snps["lineage"]
    alleles = LS_snps["Allele change"]

    # make dict like
    # {62657: ('lineage4.1', 'G/A'),
    #  107794: ('lineage4.1.2.1', 'C/T')}
    print("making positions dict...")
    pos_dict = dict(zip(snp_positions, list(zip(lins, alleles))))

    sample_dict = defaultdict(list)
    # if a file has lineage specific SNP put it to the dict (key=sample)
    print("making samples dict...")
    for file in files:
        for pos in list(snp_positions):
            with open(file) as f:
                str_to_look = str(pos) + "\t.\t" + pos_dict[pos][1][0] + "\t" + pos_dict[pos][1][2]
                if str_to_look in f.read():
                    sample_dict[file[:-7]].append(pos_dict[pos][0])

    # Filter out samples that have different lineages names (check only the first number) and
    # don't care of names "lineage4" and "lineage4.9"
    # put results in a dict
    if filter:
        print("filtering...")
        clean_and_bad = filter_mixed_samples(sample_dict, "lineage4", "lineage4.9")
        sample_dict = clean_and_bad[0]
        mixed_samples = clean_and_bad[1]
        # also you can write bad samples
        with open("mixed_samples.txt", "w") as out:
            for bad in mixed_samples:
               out.write("%s\n" % bad)

    #  And now you should reverse the sample dict
    # in such way that lineages become a keys and samples goes to lists (list are values of the dict)
    # D = {"a":["1", "2", "3"], "b":["4", "5", "3"]}
    # D_reversed = {"1":["a"], "2":["a"], "3":["a", "b"], "4":["b"], "5":["b"]}

    # reverse the sample dict
    print("reversing the sample dict...")
    lin_dict = dict_reverse(sample_dict)

    # and make from it a DataFrame of proper format
    # where the samples represents first column
    # column names are lineages
    # the "+" sign in a cell of the DataFrame indicates that this particular sample (row)
    # belongs to this particular lineage (column)
    # else NaN
    print("making a data frame...")
    new_dict = defaultdict(list)
    samples = [name[:-7] for name in files]
    for sample in samples:
        for key in lin_dict.keys():
            if sample in lin_dict[key]:
                new_dict[key].append("+")
            else:
                new_dict[key].append(np.nan)

    new_dict["samples"] = samples
    new_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in new_dict.items()]))

    # move 'samples' column to the front
    cols = new_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    new_df = new_df[cols]
    return new_df


if __name__ == '__main__':
    help_message = "At least two arguments needed:\n1 - path to phylo table, like ~/Documents/myco/compensatory/phylogenetic_snps_62.csv\n2 - output file name"
    if len(sys.argv) < 3:
        print(help_message)
        sys.exit()

    phylo_table_path = sys.argv[1] # "~/Documents/myco/compensatory/phylogenetic_snps_62.csv"
    output = sys.argv[2]
    phylo_table = make_phylo_table(ph_table=phylo_table_path, filter=False)
    phylo_table.to_csv(output, index=False)
    print("Done!")