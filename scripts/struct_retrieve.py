import os
import subprocess
import numpy as np
import pandas as pd
import argparse as ap 
from Bio.PDB import *
import seaborn as sns
from Bio import SeqIO
from pathlib import Path
from Bio.PDB import PDBList
import matplotlib.pyplot as plt


def parse_arguments(parser=None): 
    """Parse arguments given by the terminal for PDB ID."""
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args   

def plotter(in_file): 
    """
    Plot the result.csv file on seaborn barplot
     
    Args: 
        Bonds CSV file [str]
        
    Returns: 
        Seaborn barplot
    """
    
    data = pd.read_csv(in_file)

    interaction = np.array(data['interaction'])
    value = np.array(data['value'])


    sns.lineplot(x=interaction, y=value)
    plt.show()    
                                                                                                                                                                                                                                                                            
class Pedtior: 
    """"Class to edit the PDB file and then write PDB file back"""
    def __init__(self, args): 
        """Initialize the Peditor class"""
        self.args = args
        
    def struct_retrieve(self): 
        """
        Retrieve PDB structure given the PDB Id.
        
        Args
            PDB ID
        Returns 
            PDB structure file
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
        Edit the PDB file to make it ready for 
        Perl analysis. 
        Changes made: 
            1 Remove CONECT
            2 Remove HETATM HOHS
            3 Change HETTATM ATPS to LIG
        
        Args: 
            PDB file [str]
        Returns 
            Record list with changes [list]  
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
                    
    def printer(self, args, record): 
        """
        Write the record list from the 
        editor function to file. 
        
        Args: 
            Record [list]
            Address to write file [str]
        Returns: 
            PDB file written [pdb file]
        """
        with open(f"{self.args.id_input}.pdb", "w") as f: 
            for i in record: 
                f.write(i)
                
        print("\nFile written. Ready for data analysis.")
        
        
    def execute_vmd(self): 
        """
        Open the VMD tool that is present in the same directory 
        to execute the calculations for the given pdb structure 
        
        Args: 
            PDB ID [str]
            
        Returns: 
            Open VMD. 
            
            If yes: 
                Execute the VMD command script
            if no: 
                Print an error message
        """
        
        pdb_id = self.args.id_input
        
        cmd = list((f"vmd {pdb_id}.pdb source command").split(" "))
        os.chdir("vmd")
        r = subprocess.Popen(cmd)
        
        if r.communicate(): 
            return ("\nCommand executed. Now execute \
                    VMD script on terminal")
        else: 
            return ("\nPerl analysis script failed.")

args = parse_arguments()    
c = Pedtior(args)
# c.struct_retrieve()

# a = c.editor()
# c.printer(args, a)

c.execute_vmd()
#c.execute_vmd_script()

plotter("/home/nadzhou/Desktop/bonds.csv")