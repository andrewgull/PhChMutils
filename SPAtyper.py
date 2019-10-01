#!/usr/bin/env python3
# arg1: query
# arg2: min length of the repeat


import csv
import sys
import os
import glob
from datetime import datetime
import multiprocessing as mp


def find_longest_tandem(csv_tab):
    # make list of lists. Each item is a list of repeat names that follow condition:
    # diff = 1 or (diff != 1 & prev.diff = 1)
    # return the longest item (i.e. the longest list of repeat names)
    # so, first:
    # make list of differences (i.e. distances to the next repeat)
    diffs = [int(csv_tab[i + 1][8]) - int(csv_tab[i][9]) for i in range(len(csv_tab) - 1)]
    # equal its length with length of csv_tab
    diffs.append(0)
    all_tandems = list()
    bad_index = list()
    for i in range(len(csv_tab)):
        if diffs[i] == 1:
            all_tandems.append(csv_tab[i][0])
        elif diffs[i] != 1 and diffs[i - 1] == 1:
            all_tandems.append(csv_tab[i][0])
            bad_index.append(i)
        else:
            bad_index.append(i)
            all_tandems.append(csv_tab[i][0])

    if len(bad_index) == 0:
        print("found only tandem repeats")

        return all_tandems
    else:
        # split all_tandems into parts by bad_index
        all = list()
        bad_index = [-1] + bad_index

        for j in range(len(bad_index)-1):
            # if j is the last one
            if j == len(bad_index)-2:
                tandem = all_tandems[bad_index[j]+1:]
            else:
                tandem = all_tandems[bad_index[j]+1:bad_index[j+1]+1]
            all.append(tandem)

        # find the longest
        lens = list(map(lambda x: len(x), all))

        # and return the longest part
        try:
            out = [t for t in all if len(t) == max(lens)][0]
        except IndexError:
            print("no tandems found?")
            out = [""]

        return out


help_message = "Script to run blastn against all the fasta files in a directory (except the query you have to specify)" \
               "and parsing its output to get a SPA-type\n" \
               "Two arguments required:\n" \
               "1 - query name (a collection of SPA repeats in fasta format)\n" \
               "2 - min length of a repeat"

if len(sys.argv) < 3:
    print(help_message)
    sys.exit()

query = sys.argv[1]
minlen = sys.argv[2]
threads = mp.cpu_count()
filter_by_start = False

# list for output
out_list = []
# subjects list
subs = glob.glob("*.fasta")
# and remove fasta file of the query
subs = [item for item in subs if item != query]

# run blast for each subject
for subject in subs:
    print("working on %s" % subject)
    # print(subject)
    # check database exits
    db = glob.glob("%s.nhr" % subject)
    if len(db) == 0:
        print("make blast database")
        os.system("makeblastdb -in %s -dbtype nucl" % subject)
    else:
        print("blast database for %s already exists" % subject)

    # run blastn with word_size == minlen
    time = str(datetime.now())[:-7]
    time = time.split(" ")
    time = "_".join(time)
    out_name = "blastn_result_" + subject + "_" + time + ".tsv"
    command = "blastn -db %s -query %s -outfmt 6 -out %s -evalue 0.0001 -task megablast -word_size %s -num_threads %s" \
              % (subject, query, out_name, minlen, threads)
    os.system(command)

    # read blast output and sort by 9th column (i.e. start in subject)
    reader = csv.reader(open(out_name), delimiter="\t")
    # str to int in row[8]
    sorted_blast = sorted(reader, key=lambda row: int(row[8]), reverse=False)

    # value in 'q_start' column (#6 pythonic) have to be equal 1
    if filter_by_start:
        sorted_blast = [row for row in sorted_blast if row[6] == "1"]

    # find longest sequence of tandem arranged repeats
    r_types = find_longest_tandem(sorted_blast)
    # r_types = [item[0] for item in long_tandem]
    r_types = "-".join(r_types)

    output_line = subject + "\t" + r_types
    out_list.append(output_line)

# write all the results to a file
with open("analysis_result.tsv", "w") as out:
    for line in out_list:
        out.write("%s\n" % line)
print("Analysis finished successfully.\nSee 'analysis_results.tsv' for details")
# Windows CMD finish
# x = input("Done! See 'analysis_result.tsv'\nPress any key to exit")
