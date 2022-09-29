import pandas as pd 
import argparse
import matplotlib.pyplot as plt


def plot_distribution(df, study_name):
    '''
    plot the distribution of expression for the given df
    '''
    # df.plot.scatter(x="orf_count", y="cds_count", c="ratio", colormap='viridis')
    df = df.query(f'file.str.startswith("{study_name}")',  engine="python")
    print(df)
    orf_count = list(df.orf_count) 
    cds_count = list(df.cds_count)
    annotations = list(df.file)


    plt.figure(figsize=(8,6))
    plt.scatter(orf_count,cds_count)
    plt.xlabel("ORF Count")
    plt.ylabel("CDS Count")
    plt.title(f"{list(df.transcript)[0]}",fontsize=15)
    for i, label in enumerate(annotations):
        plt.annotate(label, (orf_count[i], cds_count[i]))

    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", help="path to transcript probabilites file")
    parser.add_argument("-s", help="Study first author name")
    args = parser.parse_args()
    df = pd.read_csv(args.f)
    plot_distribution(df, args.s)