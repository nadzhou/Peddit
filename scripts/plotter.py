import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

def plotter(in_file): 
    """
        Plot the bonds.csv
    """
    
    data = pd.read_csv(in_file)

    print(data.head())
    ser = np.array(data['ser'])
    interaction = np.array(data['ser1'])
    value = np.array(data['value'])
    amino_acid = np.array(data['amino_acid'])
    z = np.array([value, interaction])
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(amino_acid, value, z, cmap='binary')
    plt.show()  


def main(): 

    plotter("/home/nadzhou/DEVELOPMENT/p_editor/scripts/bonds.csv")


if __name__ == '__main__': 
    main()