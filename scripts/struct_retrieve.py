import argparse as ap 
from Bio.PDB import PDBList
from Bio.PDB import *
from pathlib import Path
from Bio import SeqIO
import subprocess

def parse_arguments(parser=None): 
    """Parse arguments given by the terminal for PDB ID."""
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args
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
        Pdb_id = self.args.id_input
        pdbl = PDBList()
        ppb = PPBuilder()

        pdbl.retrieve_pdb_file(Pdb_id, file_format='pdb', pdir=".")
        p = Path(f'pdb{Pdb_id}.ent')
        p.replace(f'{Pdb_id}.pdb')
        print("PDB file written.")

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
                if "HETATM" in line: 
                    if "HOH" not in line: 
                        record.append(line)
                    elif "ATP" in line: 
                        line = line.replace("ATP", 'LIG')
                        record.append(line)
                elif 'CONECT' not in line: 
                    record.append(line)                
        print("ATPs changed to LIG.\nHETATM water removed")
        print("File written to be written.")          
        return record
                    
    def printer(self, self.args, record): 
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
                
        print("File written. Ready for data analysis.")
                
        def execute_pl(self): 
            cmd = list(("perl analyze.pl").split(" "))
            r = subprocess.Popen(cmd)
            
            if r.communicate(): 
                return ("Command executed")
            else: 
                return ("Perl analysis script failed.")

args = parse_arguments()    
struct_seq_retrieve(args)
seq_extract(args)