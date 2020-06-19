# Credit to Viktor Demenev: https://github.com/Vikdemen
# who wrote the bondwriter Python file to convert the 
# TCL script to Python. 


# Credits: https://github.com/Eigenstate/vmd-python
# Thanks to Robin Betz for creating an installable 
# Python 3 version of VMD. 
# My tool is just a plugin on top of the VMD-Python 
# without which this Viktor's code wouldn't run. 


import vmd
from typing import List


def main():
    data = calculate_fbonds("4qdi.pdb")
    write_to_file(data, "bonds.csv")


def calculate_fbonds(filename: str) -> List[str]:
    data = []
    # noinspection PyUnresolvedReferences
    top = vmd.Molecule.Molecule().load(filename, 'pdb')
    # noinspection PyUnresolvedReferences
    resname_LIG_bonds = list(vmd.atomsel(selection="resname LIG", molid=top))
    chosen_bonds = []
    for atom_b_index in resname_LIG_bonds:
        # noinspection PyUnresolvedReferences
        close_bonds = list(vmd.atomsel(selection=f"within 5 of index {atom_b_index}", molid=top))
        for close_bond in close_bonds:
            if close_bond not in resname_LIG_bonds:
                chosen_bonds.append(close_bond)
        for atom_a_index in chosen_bonds:
            # noinspection PyUnresolvedReferences
            atom_a_selection = vmd.atomsel(selection=f'index {atom_a_index}', molid=top)
            # noinspection PyUnresolvedReferences
            distance, = vmd.measure.bond(atom_a_index, atom_b_index)
            if distance < 4:
                atom_b_selection = vmd.atomsel(selection=f"index {atom_b_index}", molid=top)
                pt1a, = atom_b_selection.name
                pt1r, = atom_b_selection.resid
                pt1rn, = atom_b_selection.resname
                pt2a, = atom_a_selection.name
                pt2r, = atom_a_selection.resid
                pt2rn, = atom_a_selection.resname

                data.append(f"{pt1r},{pt1rn},{pt1a},{pt2r},{pt2rn},{pt2a},{distance}\n")
        chosen_bonds = chosen_bonds[1:]
    return data


header = "ser,lig,pg,serl,amino_acid,interaction,value\n"


def write_to_file(data: List[str], filename):
    with open(filename, "w") as fbonds:
        fbonds.write(header)
        fbonds.writelines(data)


if __name__ == '__main__':
    main()
