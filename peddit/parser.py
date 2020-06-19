import argparse as ap

def parse_arguments(parser: 'argparse obj'=None) -> 'argparse obj': 
    """ Parse arguments given by the terminal for PDB ID.

        Returns: 
            args [argparse obj]: PDB ID argparse object
    """
    if not parser: 
        parser = ap.ArgumentParser()

    parser.add_argument("id_input", help="Input PDB ID")
    parser.add_argument("output_path", help="Output directory for PDB and bonds file.")

    args = parser.parse_args()

    return args   
