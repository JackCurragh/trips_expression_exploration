# Calculating probability of initiation and reinitiation in Trips Files 

import os
import argparse
from core import *




def calculate_studies_probabilities(read_file_paths, row, study_name=None):
    '''
    calculate uORF intiiation probabilities for a given study
    '''
    counter = 0
    probs = []
    for index, file in enumerate(read_file_paths):
        if study_name == None:
            pass
        elif study_name in file:
            pass
        else:
            continue
        read_file_path = read_file_paths[file]
        ordered_position_counts = get_unambig_reads_for_transcript(read_file_path, row['transcript'])

        if 'None' in list(row) or ordered_position_counts == None:
            print(list(row))
            continue

        cds_count = get_counts_in_range(ordered_position_counts, int(row['cds_start']), int(row['cds_stop']))
        cds_length = int(row['cds_stop']) - int(row['cds_start'])
        cds_rpb = cds_count/cds_length

        rest_of_tx_count_cds = get_counts_in_range(ordered_position_counts, int(row['cds_start']), int(row['length']))
        rest_of_tx_count_cds_length = int(row['length']) - int(row['cds_start'])
        rest_of_tx_rpb_cds = rest_of_tx_count_cds/rest_of_tx_count_cds_length

        orf_count = get_counts_in_range(ordered_position_counts, int(row['orf_start']), int(row['orf_stop']))
        orf_length = int(row['orf_stop']) - int(row['orf_start'])
        orf_rpb = orf_count/orf_length

        rest_of_tx_count_orf = get_counts_in_range(ordered_position_counts, int(row['orf_start']), int(row['length']))
        rest_of_tx_count_orf_length = int(row['length']) - int(row['orf_start'])
        rest_of_tx_rpb_orf = rest_of_tx_count_orf/rest_of_tx_count_orf_length



        if cds_count > 5 and orf_count > 5:
            ratio = round(cds_rpb/orf_rpb, 2)
            prob_orf_init = round(orf_rpb/rest_of_tx_rpb_orf, 2)
            prob_cds_init = round(cds_rpb/rest_of_tx_rpb_cds, 2)
            prob_orf_init_raw = round(orf_count/rest_of_tx_count_orf, 2)
            prob_cds_init_raw = round(cds_count/rest_of_tx_count_cds, 2)
            init_ratio = round(prob_cds_init/prob_orf_init, 2)
            probs.append([row['transcript'], file,
                                    orf_count, 
                                    cds_count, 
                                    round(ratio, 2), 
                                    prob_orf_init, 
                                    prob_cds_init, 
                                    prob_orf_init_raw, 
                                    prob_cds_init_raw, 
                                    init_ratio, 
                                    round(init_ratio - ratio,2)])
        else:
            print(row['transcript'], cds_count, orf_count)
            counter += 1
    
    
    return pd.DataFrame(probs, columns=[
                            'transcript',
                            'file',
                            'orf_count',
                            'cds_count',
                            'ratio',
                            'prob_uORF_init_rpb',
                            'prob_cds_init_rpb',
                            'prob_uORF_init_raw',
                            'prob_cds_init_raw',
                            'initiation_ratio',
                            'ratio_diff_init_minus_full'
                                ])
            


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
    single_transcript_single_orf = get_single_transcript_single_orf(openprot)

    print(single_transcript_single_orf.shape)

    print(len(riboseq_file_path))
    return False
    probs_df = None
    counter = 0 
    for row in single_transcript_single_orf.iterrows():
        print()

        if type(probs_df) == None:
            probs_df = calculate_studies_probabilities(riboseq_file_paths, row[1])  
            print(row[1]['transcript'])
      
        else:
            new_probs_df = calculate_studies_probabilities(riboseq_file_paths, row[1])
            probs_df = pd.concat([probs_df, new_probs_df])
            print(row[1]['transcript'], probs_df.shape)
        
        probs_df.to_csv("probabilites_single_tx_single_orf.csv")
    return True 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--openprot", help="path to list of genes/transcripts ")
    parser.add_argument("-s", help="path to organism sqlite ")
    parser.add_argument("-t", help="path to main trips.sqlite")

    args = parser.parse_args()
    main(args)

# res = cursor.execute(f'SELECT gene from transcripts WHERE transcript = {i}').fetchall()

# python scripts/initiaition_probability_calculations.py --openprot data/full_fiveprime_candidates.tsv -s /home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/homo_sapiens.Gencode_v25.sqlite -t /home/DATA/www/tripsviz/tripsviz/trips.sqlite