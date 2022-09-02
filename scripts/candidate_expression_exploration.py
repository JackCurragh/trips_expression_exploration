import pandas as pd 
import argparse


def plot_distribution(df):
    '''
    plot the distribution of expression for the given df
    '''

    df.plot.hist()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", help="path to transcript probabilites file")
    
    args = parser.parse_args()
    df = pd.read_csv(args.f)
    plot_distribution(df)