from ete3 import NCBITaxa
import argparse
import os
import sys
import KrakenFunctions as KF

# Initialize ete3's NCBI taxonomy object


def main():
    parser = argparse.ArgumentParser(
        description="""
        Parses Kraken2 output and check for consistency in assignment for reads. This is specially important for those reads which the assignation by kraken is 0.
        It LCA's all non-zero kmer assignations, and provides a level of consistency. 
        - Input: Kraken "classifiedreads" file
        - Output: TSV file.
        """
    )
    
    parser.add_argument("-i", help="Path to the GEMBASES history file")
    parser.add_argument("-o", help="Path to output TSV file ")
    
    args = parser.parse_args()
    ncbi = NCBITaxa()

    with open(args.i, "rt") as input, open(args.o, "w") as output:
        for gembasesID, name in input.strip.split("\t"):
            taxid,lineage,linnames = KF.extract_taxid_taxline(name, ncbi)
            output.write(f"{gembasesID}\t{name}\t{taxid}\t{lineage}\t{linnames}")

if __name__ == "__main__":
    main()    