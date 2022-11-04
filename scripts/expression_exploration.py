'''
The goal for this script is to provide expression statistics for genes/transcripts on Trips-Viz

The script works with processed reads in Sqlite files and records basic expression information on user specified transcripts 
'''


import pandas as pd
import argparse
import os 

import matplotlib.pyplot as plt

from core import *



def quantify_studys_expression(counts_df, read_file_path, filename):
    '''
    get the read counts for openprot regions for all transcripts in open prot file
    '''
    counts_df.is_copy = False
    for index, row in counts_df.iterrows():
        ordered_position_counts = get_unambig_reads_for_transcript(read_file_path, row['transcript'])
        if ordered_position_counts == None or row['cds_start'] == 'None' or row['cds_stop'] == 'None':
            continue

        counts_df.loc[index, f'{filename}_cds_count'] = get_counts_in_range(ordered_position_counts, int(row['cds_start']), int(row['cds_stop']))
        counts_df.loc[index, f'{filename}_orf_count'] = get_counts_in_range(ordered_position_counts, int(row['orf_start']), int(row['orf_stop']))
        counts_df.loc[index, f'{filename}_leader_count'] = get_counts_in_range(ordered_position_counts, int(0), int(row['cds_start']))
        counts_df.loc[index, f'{filename}_trailer_count'] = get_counts_in_range(ordered_position_counts, int(row['cds_stop']), int(row['length']))
        counts_df.loc[index, f'{filename}_transcript_count'] = get_counts_in_range(ordered_position_counts, int(0), int(row['length']))

        if counts_df.at[index, f'{filename}_orf_count'] > 0  and counts_df.at[index, f'{filename}_cds_count'] > 0:
            ratio = counts_df.at[index, f'{filename}_cds_count']/counts_df.at[index, f'{filename}_orf_count']
        elif counts_df.at[index, f'{filename}_cds_count'] > 0:
            ratio = 1
        else:
            ratio = 0

        counts_df.loc[index, f'{filename}_ratio'] = ratio
    return counts_df


def restructure_transcript_data(transcript, openprot_w_counts):
    '''
    extract a data frame for a given transcript from openprot_w_counts
    rows are files and columns relate to transcript expression
    '''
    transcript_df = openprot_w_counts[openprot_w_counts.transcript == transcript]
    data = []
    for column in transcript_df.columns:
        if str(column).endswith('cds_count'):
            filename = '_'.join(column.split('_')[:-2])

            idx = transcript_df[f'{filename}_cds_count'].index[0]

            cds_count = transcript_df[f'{filename}_cds_count'][idx]
            orf_count = transcript_df[f'{filename}_orf_count'][idx]
            ratio = transcript_df[f'{filename}_ratio'][idx]
            leader_count = transcript_df[f'{filename}_leader_count'][idx]
            trailer_count = transcript_df[f'{filename}_trailer_count'][idx]
            transcript_count = transcript_df[f'{filename}_transcript_count'][idx]

            data.append([filename, cds_count, orf_count, ratio, leader_count, trailer_count, transcript_count])

    restructured_df = pd.DataFrame(data, columns=[
                                'file',
                                'cds_count',
                                'orf_count',
                                'ratio',
                                'leader_count',
                                'trailer_count',
                                'transcript_count',
                                ])
    return restructured_df


def plot_restructured_df(restructured_df):
    '''
    produce simple plots  to explore df 
    '''
    # restructured_df.plot(x ='cds_count', y='orf_count', kind = 'hist')
    restructured_df.plot.hist()

    plt.show()


def main(args):
    if args.s != None or args.t != None:
        organism_sqlite_cursor = get_sqlite_cursor(args.s)
        trips_sqlite_cursor = get_sqlite_cursor(args.t)

        transcriptome_list = args.s.split('/')[-1].split('.')[1]
        readfile_paths = get_readfile_paths_for_organism(trips_sqlite_cursor, transcriptome_list)

        riboseq_file_paths = {}
        for study in readfile_paths:
            for file in readfile_paths[study]['riboseq']:
                riboseq_file_paths[f"{study}_{file}"] = readfile_paths[study]['riboseq'][file]
                test_file_path = readfile_paths[study]['riboseq'][file]

        riboseq_file_paths = {i:riboseq_file_paths[i] for i in list(riboseq_file_paths.keys()) }

    if os.path.exists(args.o):
        openprot_w_counts = pd.read_csv(args.o)
        global_restructured_df = None
        for transcript in openprot_w_counts.transcript:

            restructured_df = restructure_transcript_data(transcript, openprot_w_counts)
            if type(global_restructured_df) == None:
                global_restructured_df = restructured_df.assign(file = f"{transcript}_" + restructured_df.file)
            else:
                new_restructured_df = restructured_df.assign(file = f"{transcript}_" + restructured_df.file)
                global_restructured_df = pd.concat([global_restructured_df, new_restructured_df])

            # transcript_df = openprot_w_counts[openprot_w_counts.transcript == transcript][["transcript", "gene", "cds_start","cds_stop","orf_start","orf_stop"]]
            # plot_restructured_df(restructured_df)

        print(global_restructured_df.describe())
        print(global_restructured_df.sort_values(by=["orf_count"]))
    
    else:
        openprot = read_openprot_annotations(args.openprot)
        single_tx_openprot = openprot[openprot.tx_count_for_gene == 1]
        
        gene_count_df = single_tx_openprot['gene'].value_counts().to_frame()
        single_orf_genes_df = gene_count_df[gene_count_df.gene == 1]
        single_orf_genes_list = list(single_orf_genes_df.index.to_list())


        counts_df = single_tx_openprot[single_tx_openprot['gene'].isin(single_orf_genes_list)].copy()

        for idx, file in enumerate(riboseq_file_paths):
            print(file, idx, idx/len(riboseq_file_paths))
            counts_df = quantify_studys_expression(counts_df, riboseq_file_paths[file], file)
            counts_df.to_csv(args.o)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-l", help="path to list of genes/transcripts ")
    parser.add_argument("-s", help="path to organism sqlite ")
    parser.add_argument("-t", help="path to main trips.sqlite")
    parser.add_argument("--openprot", help="output path for the report")

    parser.add_argument("-o", help="output path for the report")

    args = parser.parse_args()
    main(args)


'''
python expression_exploration.py -l transcript_list.txt \
-s /home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/homo_sapiens.Gencode_v25.sqlite \
-t /home/DATA/www/tripsviz/tripsviz/trips.sqlite \
--openprot data/new_openprot_w_counts.csv \
-o example_output.txt
'''
