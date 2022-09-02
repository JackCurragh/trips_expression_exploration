import numpy as np
import argparse
from core import *
from weblogo import *
import pandas as pd
import matplotlib.pyplot as plt


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
    probs_df.reset_index(drop=True, inplace=True)

    tx_grouped = probs_df.groupby(by="transcript") 
    tx_grouped = tx_grouped.mean()

    init = []
    for row in tx_grouped.iterrows():

        openprot_subset = openprot[openprot.transcript == row[0]]
        idx = openprot_subset.index[0]
        orf_start = openprot_subset['orf_start'][idx]
        orf_stop = openprot_subset['orf_stop'][idx]

        cds_start = openprot_subset['cds_start'][idx]
        cds_stop = openprot_subset['cds_stop'][idx]

        init.append([f'{row[0]}_orf', orf_start, orf_stop, row[1]["prob_uORF_init_raw"]])
        init.append([f'{row[0]}_cds', cds_start, cds_stop, row[1]["prob_cds_init_raw"]])

    init_df = pd.DataFrame(init, columns=[
                                'transcript',
                                'start',
                                'stop',
                                'raw_prob',
    ])
    bins = np.linspace(init_df.raw_prob.min(), init_df.raw_prob.max(), 10)
    groups = init_df.groupby(np.digitize(init_df.raw_prob, bins))

    print(init_df.shape)
    split_dfs = np.array_split(init_df, 10)

    for df in split_dfs:
        print(df.shape)
    # for group in groups:
    #     print(group[1])
    # print(groups.mean())
    # print(groups.median())

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