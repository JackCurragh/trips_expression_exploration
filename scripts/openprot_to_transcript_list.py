import pandas as pd 
import sqlite3
import argparse

def openprot_to_df(path_to_tsv, n_row_limit=1000000000000):
    '''
    read in teh downloaded openprot file and return a df wiht correct column names 
    '''

    df = pd.read_csv(path_to_tsv, sep='\t', skiprows=1, nrows=n_row_limit)
    return df


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



def main(path_to_tsv, path_to_sqlite, outfile, localization):
    openprot_df = openprot_to_df(path_to_tsv)
    transcript_df = openprot_df[['transcript accession', 'start transcript coordinates', 'stop transcript coordinates', 'localization']]

    num_rows = len(transcript_df.index)

    with open(outfile, 'w') as f:
        f.write(f"transcript\tgene\tlength\tcds_start\tcds_stop\torf_start\torf_stop\ttx_count_for_gene\n")
        sqlite_cursor = get_sqlite_cursor(path_to_sqlite)
        for index, row in transcript_df.iterrows():
            if index % 100 == 0:
                print(index, round(index/num_rows, 2))
            if row['localization'] == localization:
                transcipt_name, transcript_version = str(row['transcript accession']).split('.')[0], str(row['transcript accession']).split('.')[1]
                information_dict = get_transcript_annotation(sqlite_cursor, transcipt_name)
                information_dict['orf_start'] = int(row['start transcript coordinates'])
                information_dict['orf_stop'] = int(row['stop transcript coordinates'])
                count = get_transcript_count_for_gene(information_dict['gene'], sqlite_cursor)
                
                f.write(f"{information_dict['transcript']}\t{information_dict['gene']}\t{information_dict['length']}\t{information_dict['cds_start']}\t{information_dict['cds_stop']}\t{information_dict['orf_start']}\t{information_dict['orf_stop']}\t{count}\n")

            else: 
                continue        



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", help="path to organism sqlite ")
    parser.add_argument("-t", help="path to openprot tsv")
    parser.add_argument("-o", help="output path for the report")

    args = parser.parse_args()
    main(args.t, args.s, args.o, "5'UTR")


    '''

python openprot_to_transcript_list.py -s /home/DATA/www/tripsviz/tripsviz/trips_annotations/homo_sapiens/homo_sapiens.Gencode_v25.sqlite -t /home/gwips/projects/trips_expression_exploration/human-openprot-r1_6-altprots+isoforms-grch38.95.tsv -o 5_prime_candidates.tsv
    '''