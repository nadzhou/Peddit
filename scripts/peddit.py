import os
import subprocess
import numpy as np
import pandas as pd
import seaborn as sns

from Bio import SeqIO
import argparse as ap 
from Bio.PDB import *
from pathlib import Path

from Bio.PDB import PDBList
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def parse_arguments(parser=None): 
    """Parse arguments given by the terminal for PDB ID.
    """
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args   


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
                                                                                                                                                                                                                                                                            
class Pedtior: 
    """"Class to edit the PDB file and then write PDB file back
    """

    def __init__(self, args): 
        """Initialize the Peditor class"""
        self.args = args
        
    def struct_retrieve(self): 
        """
            Retrieve PDB structure given argparse ID
        """
        Pdb_id = str(self.args.id_input)
        pdbl = PDBList()
        ppb = PPBuilder()

        retrieve = pdbl.retrieve_pdb_file(Pdb_id, file_format='pdb', pdir=".")
        if retrieve:  
            p = Path(f'pdb{Pdb_id}.ent')
            p.replace(f'{Pdb_id}.pdb')

            print("\nPDB file written.")

    def editor(self): 
        """
            Edit the PDB file. Remove HOH, change ATP to LIG. 
        """
        record = []
        with open(f"{self.args.id_input}.pdb", "r") as f: 
            for line in f: 
                if 'HETATM' in line: 
                    if "ATP" in line: 
                        line = line.replace('ATP', 'LIG')
                        record.append(line)
                    if 'HOH' not in line: 
                        record.append(line)
                    
                elif 'CONECT' not in line: 
                    record.append(line)                
        print("ATPs changed to LIG.\nHETATM water removed")
        print(f"\n Edited file has been rewritten to {self.args.id_input}")          
        return record


    def edited_pdb_writer(self, args, record): 
        """
            Write the edited PDB list to file. 
        """
        with open(f"{self.args.id_input}.pdb", "w") as f: 
            for i in record: 
                f.write(i)
                
        print("\nFile written. Ready for data analysis.")
        
        
    def execute_vmd(self): 
        """
            Open up VMD and execute command. 
        """
        
        pdb_id = self.args.id_input
        
        cmd = list((f"vmd {pdb_id}.pdb -e command").split(" "))
        os.chdir("vmd")
        r = subprocess.Popen(cmd)
        
        if r.communicate(): 
            print("\nCommand executed. Now execute \
                    VMD script on terminal")


def main(): 
    print("Initiating Peddit tool.")
    args = parse_arguments()    

    c = Pedtior(args)
    c.struct_retrieve()

    a = c.editor()
    c.edited_pdb_writer(args, a)

    print("Executing VMD script..\n")
    c.execute_vmd()

    # Plot the output file from VMD
    plotter("/home/nadzhou/Desktop/bonds.csv")



if __name__ == '__main__': 
    main()