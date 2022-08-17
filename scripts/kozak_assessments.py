import sqlite3 
import argparse
from core import *
from weblogo import *
import pandas as pd


def get_transcript_sequence(transcript, organism_sqlite_cursor):
    '''
    return the nucleotide sequence for a given transcript as a string
    '''
    sequence = organism_sqlite_cursor.execute(
        f"SELECT sequence FROM transcripts WHERE transcript == '{transcript}'"
    ).fetchall()[0]

    return sequence[0]


def main(args):
    cursor = get_sqlite_cursor(args.s)

    openprot = read_openprot_annotations(args.transcripts)

    probs_df = pd.read_csv(args.p)
    init = []
    for row in probs_df.iterrows():
        init.append([f'{row[1]["transcript"]}_orf', row[1]["prob_uORF_init_rpb"], row[1]["prob_uORF_init_raw"]])
        init.append([f'{row[1]["transcript"]}_cds', row[1]["prob_cds_init_rpb"], row[1]["prob_cds_init_raw"]])

    init_df = pd.DataFrame(init, columns=[
                                'transcript',
                                'rpb_prob',
                                'raw_prob',
    ])
    print(init_df.sort_values('raw_prob'))
    unique_transcripts = list(probs_df['transcript'].unique())
    return True

    outfile = open(args.o, 'w')
    for row in openprot.iterrows():
        if "None" in list(row[1]):
            continue

        cds_start = int(row[1]['cds_start'])
        orf_start = int(row[1]['orf_start'])
        sequence = get_transcript_sequence(row[1]['transcript'], cursor)
        outfile.write(f">{row[1]['transcript']}_orf\n")
        outfile.write(f"{sequence[orf_start-10:orf_start+10]}\n")
        outfile.write(f">{row[1]['transcript']}_cds\n")
        outfile.write(f"{sequence[cds_start-10:cds_start+10]}\n")
    outfile.close
    return True 



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--transcripts", help="path to subset of openprot data ")
    parser.add_argument("-s", help="path to organism sqlite ")
    parser.add_argument("-p", help="path to probabilities csv")
    parser.add_argument("-o", help="path to output fasta file of sequences at start site")

    args = parser.parse_args()
    main(args)