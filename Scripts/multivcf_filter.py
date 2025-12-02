#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import argparse
import sys
import numpy as np

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

def read_multifasta(multifasta): 
    """Read multifasta file into pandas DataFrame object."""
    f = open(multifasta, "r")
    all = f.read().splitlines()
    names= [item.strip(">") for item in (all[::2])]
    seqs = (all[1::2])
    data = [list(seq) for seq in seqs]
    df = pd.DataFrame(data, index = names)
    f.close()
    return df

def percent_deletion(geno, perc_deletion):
    """Filter SNP positions based on partial deletion value (e.g. 0.95 for 95% partial deletion)."""
    snp_ct = geno.shape[0]
    #Good SNPs/ total length greater than cutoff
    bool_array = geno[geno == 'N'].isna().sum()/snp_ct >= float(perc_deletion)
    subset_df = geno.loc[:, bool_array]
    return subset_df

def check_biallelic(array):
    """Keep only biallelic SNPs."""
    tmp = array.replace('N', np.nan)
    bool_array = tmp.nunique(dropna=True) == 2
    subset_df = array.loc[:, bool_array]
    return subset_df

def write_multifasta(subset_df, output_file):
    """Write filtered output in multifasta format."""
    with open(output_file, 'w') as f: 
        for index, row in subset_df.iterrows(): 
            f.write(">" + index)
            f.write("\n")
            f.write("".join(row))
            f.write("\n")
    f.close()

def main(): 
    parser = argparse.ArgumentParser(description = "Parameters for multifasta_subset.py.")
    parser.add_argument('multifasta', help = "Input multifasta file.")
    parser.add_argument('perc_deletion', type=float, choices=[Range(0.0, 1.0)], help = "Percent deletion for filtering multifasta (e.g. 0.95 for 95%% partial deletion).")
    parser.add_argument('output_dir', help = "Directory for output files.")
    parser.add_argument('output', help = "Base prefix for naming output files.")
    try: 
        args = parser.parse_args()
    except: 
        parser.print_help()
        sys.exit(0)

    df = read_multifasta(args.multifasta)
    print ("Multifasta data loaded.")
    pd_df = percent_deletion(df, args.perc_deletion)
    biallelic_pd_df = check_biallelic(pd_df)
    print ("Partial deletion filtering complete. Outputting fasta file.")
    write_multifasta(biallelic_pd_df, args.output_dir + '/' + args.output + '.fasta')
    print ("Done!")

if __name__ == "__main__": 
    main()
