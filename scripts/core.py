#Core functions used to explore expression in Trips Files 

import sqlite3 
from sqlitedict import SqliteDict
import pickle5
import collections
import pandas as pd
import os 

def get_sqlite_cursor(db_path):
    '''
    open a connection to the organism sqlite db and cursor
    '''
    print(f'Opening database file at {db_path}')
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
    if os.path.exists(filepath):
        sqlite_db = SqliteDict(f"{filepath}", autocommit=False, decode=my_decoder)
        return sqlite_db
    
    else:
        return None 


def get_unambig_reads_for_transcript(file_path, transcipt_name):
    '''
    return the read dict for unambiguous read mappings on 'transcipt'
    '''
    position_counts = {}
    sqlite_dict = read_sqlite_dict(file_path)
    if sqlite_dict == None:
        return None

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



def get_transcript_annotation(organism_cursor, transcript):
    '''
    return transcript information from sqlite db
    '''
    try:
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
    except:
        information_dict = {
            'transcript':transcript,
            'gene':'None',
            'length':'None',
            'cds_start':'None',
            'cds_stop':'None',
        }
    return information_dict


def get_transcript_count_for_gene(gene, sqlite_cursor):
    '''
    return an integer of the number of transcripts a gene has in the db
    '''
    count = sqlite_cursor.execute(
        f"SELECT COUNT(transcript) FROM transcripts WHERE gene == '{gene}'"
    ).fetchall()[0][0]
    return count


def read_openprot_annotations(path_to_openprot_file):
    '''
    read in openprot transcript list and return a df
    '''
    df = pd.read_csv(path_to_openprot_file, sep='\t')
    df = df.drop_duplicates()
    return df


def get_single_transcript_single_orf(openprot):
    '''
    from the openprot input return a dataframe of entries with one transcript and one orf 
    '''
    single_tx_openprot = openprot[openprot.tx_count_for_gene == 1]
    
    gene_count_df = single_tx_openprot['gene'].value_counts().to_frame()
    single_orf_genes_df = gene_count_df[gene_count_df.gene == 1]
    single_orf_genes_list = list(single_orf_genes_df.index.to_list())
    single_gene_single_orf = single_tx_openprot[single_tx_openprot['gene'].isin(single_orf_genes_list)].copy()
    return single_gene_single_orf


def get_single_transcript_N_orf(openprot, N):
    '''
    from the openprot input return a dataframe of entries with one transcript and one orf 
    '''
    single_tx_openprot = openprot[openprot.tx_count_for_gene == 1]
    
    gene_count_df = single_tx_openprot['gene'].value_counts().to_frame()
    single_orf_genes_df = gene_count_df[gene_count_df.gene <= N]
    single_orf_genes_list = list(single_orf_genes_df.index.to_list())
    single_gene_single_orf = single_tx_openprot[single_tx_openprot['gene'].isin(single_orf_genes_list)].copy()
    return single_gene_single_orf