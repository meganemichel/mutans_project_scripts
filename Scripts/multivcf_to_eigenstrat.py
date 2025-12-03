#!/usr/bin/env python

from numpy import genfromtxt
import snp_data_processing
import argparse
import sys


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

def main(): 
    parser = argparse.ArgumentParser(description = "Parameters for multifasta_subset.py.")
    parser.add_argument('snpAlignment', help = "Input snpAlignment from multivcfAnalyzer.")
    parser.add_argument('snpTable', help = "Input snpTable from multivcfAnalyzer.")
    parser.add_argument('perc_deletion', type=float, choices=[Range(0.0, 1.0)], help = "Percent deletion for filtering multifasta (e.g. 0.95 for 95%% partial deletion).")
    parser.add_argument('output_dir', help = "Directory for output files.")
    parser.add_argument('output', help = "Base prefix for naming output files.")
    try: 
        args = parser.parse_args()
    except: 
        parser.print_help()
        sys.exit(0)

    print ('Reading SNP Alignment and tsv table.')
    alignment = snp_data_processing.read_multifasta(args.snpAlignment)
    snp_table = genfromtxt(args.snpTable, delimiter='\t',  dtype=str)
    snp_dist, ref_calls, genotypes = snp_data_processing.process_snp_table(snp_table)
    pd_df = snp_data_processing.percent_deletion(alignment, args.perc_deletion)
    biallelic_pd_df = snp_data_processing.check_biallelic(pd_df)
    genotypes_subset, ref_call_subset, snp_position_subset = snp_data_processing.subset_snp_table(biallelic_pd_df, snp_dist, ref_calls, genotypes)
    alt_allele_subset = snp_data_processing.get_alt_allele(genotypes_subset)
    eigen_encoding = snp_data_processing.get_geno(genotypes_subset)

    with open(args.output_dir + "/" + args.output + ".geno", 'w') as f:
        for snp in eigen_encoding:
            for item in snp:
                f.write((item))
            f.write("\n")
    f.close()

    with open(args.output_dir + "/" + args.output + ".ind", 'w') as f1:
        for name in biallelic_pd_df.index:
            f1.write(name)
            f1.write("\t")
            f1.write("U")
            f1.write("\t")
            f1.write("Unknown")
            f1.write("\n")
    f1.close()

    with open(args.output_dir + "/" + args.output + ".snp", 'w') as f2:
        for index, snp in enumerate(snp_position_subset):
            f2.write(str(int(index) + 1))
            f2.write("\t")
            f2.write("1")
            f2.write("\t")
            f2.write("0")
            f2.write("\t")
            f2.write(str(snp))
            f2.write("\t")
            f2.write(str(ref_call_subset[index]))
            f2.write("\t")
            f2.write(alt_allele_subset[index])
            f2.write("\n")
    f2.close()

    with open(args.output_dir + "/" + args.output + '.fasta', 'w') as f3: 
        for index, row in biallelic_pd_df.iterrows(): 
            f3.write(">" + index)
            f3.write("\n")
            f3.write("".join(row))
            f3.write("\n")
    f3.close()

if __name__ == "__main__": 
    main()
