import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def plotter(in_file): 
    """
        Plot the bonds.csv
    """
    
    data = pd.read_csv(in_file)
    
    fig = plt.figure()

    ax = sns.lineplot(x='amino_acid', 
                        y="interaction", 
                        data=data)

    return ax 


def main(): 

    plotter("/home/nadzhou/DEVELOPMENT/p_editor/scripts/bonds.csv")
    plt.show()

if __name__ == '__main__': 
    main()