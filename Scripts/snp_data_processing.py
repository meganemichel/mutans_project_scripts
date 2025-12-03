#!/usr/bin/env python
import pandas as pd
import numpy as np

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

def process_snp_table(snpTable):
    """Usage: Get SNP postitions, reference alleles, and genotypes for all variants in output of multivcf analyzer.
    Input: 1)array- output of multivcf analyzer read into array
    """
    snp_position = np.ndarray.flatten( snpTable[1:, [0]]).astype(int)
    ref_call = np.ndarray.flatten( snpTable[1:, [1]]).astype(str)
    genotypes = snpTable[1:, 2:]
    return snp_position, ref_call, genotypes

def subset_snp_table(subset_alignment, snp_position, ref_call, genotypes): 
    indices = list(subset_alignment.columns)
    genotypes_subset = genotypes[indices,:]
    ref_call_subset = ref_call[indices]
    snp_position_subset = snp_position [indices] 
    return genotypes_subset, ref_call_subset, snp_position_subset

def get_alt_allele(genotypes):
    z = np.array(['.','N'])
    alt_allels = []
    for index, row in enumerate(genotypes):
        a = list(np.setdiff1d(row,z))
        alt_allels.append(",".join(a))
    return alt_allels

def get_geno(array):
    array[array == "A"] = 0
    array[array == "C"] = 0
    array[array == "G"] = 0
    array[array == "T"] = 0
    array[array == "."] = 2
    array[array == "N"] = 9
    return array