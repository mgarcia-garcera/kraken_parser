import sys
from ete3 import NCBITaxa
from collections import defaultdict
import KrakenFunctions as KF
import argparse
import os
ncbi = NCBITaxa()

def main():
    parser = argparse.ArgumentParser(
        description="""
        Parses Kraken2 output and check for consistency in assignment for reads. This is specially important for those reads which the assignation by kraken is 0.
        It LCA's all non-zero kmer assignations, and provides a level of consistency. 
        - Input: Kraken "classifiedreads" file
        - Output: TSV file.
        """
    )
    
    parser.add_argument("-i", help="Path to the kraken2 \"classifiedreads\" file ")
    parser.add_argument("-o", help="Path to output TSV file ")
    
    args = parser.parse_args()

    if not args.i:
        script = os.path.basename(sys.argv[0])
        print(f"Input file Unknown or not provided\nUsage: python {script} -i <kraken2_output> -o <output.tsv>")
        sys.exit(1)
    elif not args.o:
        script = os.path.basename(sys.argv[0])
        print(f"Input file Unknown or not provided\nUsage: python {script} -i <kraken2_output> -o <output.tsv>")
        sys.exit(1)

    KF.process_kraken2(args.i, args.o,ncbi)