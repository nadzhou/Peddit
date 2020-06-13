import os
import subprocess
from Bio import SeqIO
import argparse as ap 

from pathlib import Path
from Bio.PDB import PDBList

from peddit.bondwriter import calculate_fbonds
from peddit.bondwriter import write_to_file

from pathlib import Path
from typing import List

                                                                                                                                                                                                                                                                            
class Peditor: 
    """"Class to edit the PDB file and then write PDB file back
    """

    def __init__(self, args): 
        """Initialize the Peditor class
        """
        self.args = args
        self.pdb_id = self.args.id_input
        self.out_dir = self.args.output_path
        self.out_dir = self.make_output_dir()
        
    
    def make_output_dir(self) -> 'pathlib obj': 
        self.out_dir = Path(self.out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)

        return self.out_dir

    def struct_retrieve(self): 
        """
            Retrieve PDB structure given argparse ID
        """
        pdbl = PDBList()

        pdbl.retrieve_pdb_file(self.pdb_id, file_format='pdb', pdir=f"{self.out_dir}/")
        

    def replace_ent_to_pdb_name(self): 
        """Replace the pdb-.ent name with self.pdb_id.pdb
        """
        p = Path(f'{self.out_dir}/pdb{self.pdb_id}.ent')
        p.replace(f'{self.out_dir}/{self.pdb_id}.pdb')


    def editor(self) -> List: 
        """
            Edit the PDB file. Remove HOH, change ATP to LIG. 
        """
        record = []
        with open(f"{self.out_dir}/{self.pdb_id}.pdb", "r") as f: 
            for line in f: 
                if 'HETATM' in line: 
                    if "ATP" in line: 
                        line = line.replace('ATP', 'LIG')
                        record.append(line)

                    if 'HOH' not in line: 
                        record.append(line)
                    
                elif 'CONECT' not in line: 
                    record.append(line) 

        print("\nATPs changed to LIG.\nHETATM water removed")
        print(f"\nEdited file has been rewritten to {self.args.id_input}")  

        return record


    def edited_pdb_writer(self, args: 'argparse obj', record: List[str]): 
        """
            Write the edited PDB list to file. 
        """
        with open(f"{self.out_dir}/{self.pdb_id}.pdb", "w") as f: 
            for i in record: 
                f.write(i)
                
        print("\nFile written. Ready for data analysis.")
        
        
    def calculate_interaction(self, bonds_out_file: str): 
        """ Invoke the bondwriter module to calculate bonds and write to file

            Args: 
                bonds_out_file [str]: Output CSV file

        """
        bonds_data = calculate_fbonds(f"{self.out_dir}/{self.pdb_id}.pdb")

        write_to_file(bonds_data, self.out_dir/bonds_out_file)


def main(): 
    print("Initiating Peddit tool.")
    args = parse_arguments()    

    c = Peditor(args)
    c.struct_retrieve()
    c.replace_ent_to_pdb_name()

    a = c.editor()
    c.edited_pdb_writer(args, a)
    print("\nPDB file written.")


    print(f"Calculating bonds for {args.id_input}.pdb...")
    c.calculate_interaction(f"{args.id_input}_bonds.csv")
    print("Writing bonds interaction file: done")


    print("\nPlotting the data...")
    plotter(f"{args.id_input}_bonds.csv")
    plt.show()


if __name__ == '__main__': 
    main()