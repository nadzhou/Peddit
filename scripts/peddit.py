import os
import subprocess
from Bio import SeqIO
import argparse as ap 

from pathlib import Path
from Bio.PDB import PDBList
<<<<<<< HEAD
=======

from typing import List
>>>>>>> e15c505e56d7f73deaadb53b7b7f015aee12ed43

from typing import List

def parse_arguments(parser: 'argparse obj'=None) -> 'argparse obj': 
<<<<<<< HEAD
    """  Parse arguments given by the terminal for PDB ID.

        Returns: 
            args [argparse obj]: PDB ID argparse object
=======
    """Parse arguments given by the terminal for PDB ID.
>>>>>>> e15c505e56d7f73deaadb53b7b7f015aee12ed43
    """
    if not parser: 
        parser = ap.ArgumentParser()
    parser.add_argument("id_input", help="Input PDB ID")
    args = parser.parse_args()

    return args   

                                                                                                                                                                                                                                                                            
class Pedtior: 
    """"Class to edit the PDB file and then write PDB file back
    """

    def __init__(self, args): 
        """Initialize the Peditor class
        """
        self.args = args
        

    def struct_retrieve(self): 
        """
            Retrieve PDB structure given argparse ID
        """
        self.pdb_id = str(self.args.id_input)
        pdbl = PDBList()

        pdbl.retrieve_pdb_file(self.pdb_id, file_format='pdb', pdir=".")
        

    def replace_ent_to_pdb_name(self): 
        """Replace the pdb-.ent name with self.pdb_id.pdb
        """
        p = Path(f'pdb{self.pdb_id}.ent')
        p.replace(f'{self.pdb_id}.pdb')

        print("\nPDB file written.")


<<<<<<< HEAD
=======

>>>>>>> e15c505e56d7f73deaadb53b7b7f015aee12ed43
    def editor(self) -> List: 
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


    def edited_pdb_writer(self, args: 'argparse obj', record: List[str]): 
        """
            Write the edited PDB list to file. 
        """
        with open(f"{self.pdb_id}.pdb", "w") as f: 
            for i in record: 
                f.write(i)
                
        print("\nFile written. Ready for data analysis.")
        
        
    def execute_vmd(self): 
        """
            Open up VMD and execute command. 
        """
<<<<<<< HEAD
                
        # start VMD.lnk 4xuh.pdb -e command
        cmd = list((f"vmd {self.pdb_id} -e command.txt").split(" "))
=======
        
        pdb_id = self.args.id_input
        
        # start VMD.lnk 4xuh.pdb -e command
        cmd = list((f"vmd {pdb_id} -e command").split(" "))
>>>>>>> e15c505e56d7f73deaadb53b7b7f015aee12ed43
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
    c.replace_ent_to_pdb_name()

    print("Executing VMD script..\n")
    c.execute_vmd()


if __name__ == '__main__': 
    main()