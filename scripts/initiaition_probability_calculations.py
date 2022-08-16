# Calculating probability of initiation and reinitiation in Trips Files 

import os
import argparse
from core import *



def main(args):
    '''
    run script with inputted args
    '''    
    trips_sqlite_cursor = get_sqlite_cursor(args.t)

    transcriptome_list = args.s.split('/')[-1].split('.')[1]
    readfile_paths = get_readfile_paths_for_organism(trips_sqlite_cursor, transcriptome_list)

    riboseq_file_paths = {}
    for study in readfile_paths:
        for file in readfile_paths[study]['riboseq']:
            riboseq_file_paths[f"{study}_{file}"] = readfile_paths[study]['riboseq'][file]

    riboseq_file_path = {i:riboseq_file_paths[i] for i in list(riboseq_file_paths.keys()) }

    openprot = read_openprot_annotations(args.openprot)
    single_tx_openprot = openprot[openprot.tx_count_for_gene == 1]
    
    gene_count_df = single_tx_openprot['gene'].value_counts().to_frame()
    single_orf_genes_df = gene_count_df[gene_count_df.gene == 1]
    single_orf_genes_list = list(single_orf_genes_df.index.to_list())
    single_gene_single_orf = single_tx_openprot[single_tx_openprot['gene'].isin(single_orf_genes_list)].copy()
    print(single_gene_single_orf)
    return True 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--openprot", help="path to list of genes/transcripts ")
    parser.add_argument("-s", help="path to organism sqlite ")
    parser.add_argument("-t", help="path to main trips.sqlite")

    args = parser.parse_args()
    main(args)