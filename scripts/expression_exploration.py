'''
The goal for this script is to provide expression statistics for genes/transcripts on Trips-Viz

The script works with processed reads in Sqlite files and records basic expression information on user specified transcripts 
'''

import sqlite3 
from sqlitedict import SqliteDict
import pandas as pd
import argparse
import os 
import pickle5
import collections
import matplotlib.pyplot as plt



def get_sqlite_cursor(db_path):
    '''
    open a connection to the organism sqlite db and cursor
    '''
    connection = sqlite3.connect(f"{db_path}")
    cursor = connection.cursor()
    return cursor


def get_transcript_annotation(organism_cursor, transcript):
    '''
    return transcript information from sqlite db
    '''
    reference_details = organism_cursor.execute(
        f"SELECT transcript, gene, length, cds_start, cds_stop FROM transcripts WHERE transcript == '{transcript}'"
    ).fetchall()[0]
    information_dict = {
        'transcript':reference_details[0],
        'gene':reference_details[1],
        'length':reference_details[2],
        'cds_start':reference_details[3],
        'cds_stop':reference_details[4],
    }

    return information_dict


def get_readfile_paths_for_organism(trips_sqlite_cursor, transcriptome_list):
    '''
    return the filepaths for all read files for given organism 
    '''
    readfile_paths = {}
    organism = trips_sqlite_cursor.execute(
        f"SELECT organism_name FROM organisms WHERE transcriptome_list == '{transcriptome_list}'"
    ).fetchall()[0][0]

    studies = trips_sqlite_cursor.execute(
        f"SELECT study_id, study_name FROM studies WHERE organism_id == (SELECT organism_id FROM organisms WHERE transcriptome_list == '{transcriptome_list}') AND owner = 1 OR owner = 35"
    ).fetchall()
    study_list = [(i[0], i[1]) for i in studies]

    for study in study_list:
        if study[1] not in readfile_paths:
            readfile_paths[study[1]] = {'riboseq':{}, 'rnaseq':{}}
        
        files = trips_sqlite_cursor.execute(
            f"SELECT file_type,file_name,study_id,file_id,owner FROM files WHERE study_id == '{study[0]}'"
        ).fetchall()
        for file in files:
            information_dict = {
            'file_type':file[0],
            'file_name':file[1].replace(".shelf","").replace(".sqlite",""),
            'study_id':file[2],
            'file_id':file[3],
            'owner':file[4],
            }
            path = f"/home/DATA/www/tripsviz/tripsviz/trips_shelves/{information_dict['file_type']}/{organism}/{study[1]}/{information_dict['file_name']}.sqlite"
            if information_dict['file_type'] in readfile_paths[study[1]]:
                readfile_paths[study[1]][information_dict['file_type']][information_dict['file_name']] = path
    return readfile_paths


def my_decoder(obj):
	return pickle5.loads(obj)


def read_sqlite_dict(filepath):
    '''
    read in sqlite dict from filepath
    '''
    sqlite_db = SqliteDict(f"{filepath}", autocommit=False, decode=my_decoder)
    return sqlite_db


def get_unambig_reads_for_transcript(file_path, transcipt_name):
    '''
    return the read dict for unambiguous read mappings on 'transcipt'
    '''
    position_counts = {}
    sqlite_dict = read_sqlite_dict(file_path)
    if transcipt_name in sqlite_dict:
        read_dict = sqlite_dict[transcipt_name]['unambig']
    else:
        read_dict = {}
    for read_length in read_dict:
        for position in read_dict[read_length]:
            if position not in position_counts:
                position_counts[position] = read_dict[read_length][position]
            
            else:
                position_counts[position] += read_dict[read_length][position]


    ordered_position_counts = collections.OrderedDict(sorted(position_counts.items()))
    return ordered_position_counts
    

def get_counts_in_range(ordered_position_counts, start, stop):
    '''
    return the total counts mapped to a range 
    '''
    total_count = 0

    for position, count in ordered_position_counts.items():
        if position >= start and position <= stop:
            total_count += count

    return total_count


def read_openprot_annotations(path_to_openprot_file):
    '''
    read in openprot transcript list and return a df
    '''
    df = pd.read_csv(path_to_openprot_file, sep='\t')
    df = df.drop_duplicates()
    return df


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
            cds_count = transcript_df[f'{filename}_cds_count'][0]
            orf_count = transcript_df[f'{filename}_orf_count'][0]
            ratio = transcript_df[f'{filename}_ratio'][0]
            data.append([filename, cds_count, orf_count, ratio])

    restructured_df = pd.DataFrame(data, columns=[
                                'file',
                                'cds_count',
                                'orf_count',
                                'ratio'
                                ])
    return restructured_df


def plot_restructured_df(restructured_df):
    '''
    produce simple plots  to explore df 
    '''
    # restructured_df.plot(x ='cds_count', y='orf_count', kind = 'hist')
    restructured_df.plot.kde()

    # plt.show()


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

    # for transcript in openprot_w_counts.transcript:
    #     restructured_df = restructure_transcript_data(transcript, openprot_w_counts)
    #     plot_restructured_df(restructured_df)
    #     break

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
python expression_exploration.py -l transcript_list.txt -s /home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/homo_sapiens.Gencode_v25.sqlite -t /home/DATA/www/tripsviz/tripsviz/trips.sqlite -o example_output.txt
'''
