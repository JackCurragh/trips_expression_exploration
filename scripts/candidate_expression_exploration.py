import pandas as pd 
import argparse
import matplotlib.pyplot as plt


def plot_distribution(df, study_name):
    '''
    plot the distribution of expression for the given df
    '''
    # df.plot.scatter(x="orf_count", y="cds_count", c="ratio", colormap='viridis')
    # df = df.query(f'file.str.startswith("{study_name}") or file.str.startswith("Calv")',  engine="python")
    # df = df.query(f'not file.str.endswith("310") and not file.str.endswith("312") and file.str.startswith("Lintner17")',  engine="python")

    orf_count = list(df.orf_count) 
    cds_count = list(df.cds_count)
    annotations = list(df.file)
    print(df)
    plt.figure(figsize=(8,6))
    plt.scatter(orf_count,cds_count)
    plt.xlabel("ORF Count")
    plt.ylabel("CDS Count")
    plt.title(f"{list(df.transcript)[0]}",fontsize=15)
    for i, label in enumerate(annotations):
        plt.annotate(label, (orf_count[i], cds_count[i]))

    plt.show()

def explore_df(df):
    study_dict = {}

    for index, row in df.iterrows():
        if row.file.split('_')[0] not in study_dict:
            study_dict[row.file.split('_')[0]] = [(row.orf_count, row.cds_count, row.file.split('_')[1])]
        else:
            study_dict[row.file.split('_')[0]].append((row.orf_count, row.cds_count, row.file.split('_')[1]))

    res = []
    for i in study_dict:
        if len(study_dict[i]) < 2:
            continue
        diff_counts = {i[2]:i[1] - i[0] for i in study_dict[i]}
        if min(diff_counts.values()) < 0 and max(diff_counts.values())>100: 
            print(diff_counts)
            plot_distribution(df, i)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", help="path to transcript probabilites file")
    
    args = parser.parse_args()
    df = pd.read_csv(args.f)
    # explore_df(df)
    plot_distribution(df, "Park16")