from peddit.parser import parse_arguments
from peddit.peddit import Peditor
from peddit.plotter import plotter

import matplotlib.pyplot as plt


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
    plotter(f"{args.output_path}/{args.id_input}_bonds.csv")
    plt.show()


if __name__ == '__main__': 
    main()