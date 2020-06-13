# Peddit - A PDB editor plugin

_New_: Added vmd-python ([Robin Betz](https://github.com/Eigenstate)) and bondwriter ([Viktor](https://github.com/Vikdemen)) 
to make the tool indepdendent of the binary VMD to installable vmd-python tool. 

* Now the code will run as a Python script instead of a VMD script. 

Python script for PDB editing and analysis. Tool is a handy addition to automatically 
retrieveing PDB structures which is then used to make calculations on via VMD. 

![image](https://user-images.githubusercontent.com/25282805/77248658-f287d280-6c5c-11ea-8e22-9b1dbc140b38.png)

## MOTIVATION
It is hard to manually look through a PDB file and change sets of strings manually and 
also to then port the PDB file to VMD, the software. This tool fixes that issue, and then 
helps visualize the data that has been output on seaborn barplot. 
Thus automating the whole process, this tool adds to my arsenal of piping dev that 
I've been doing for quite a while now. 

![image](https://user-images.githubusercontent.com/25282805/77248689-1814dc00-6c5d-11ea-85fa-94f4321de1d4.png)

## USAGE
On the command line: 

    $ python main.py <pdb_idb> <output directory> 

This will run the code and plot the data shown below as a result. 

Example: 

    $ python main.py 1fmw /home/1fmw

## DEPENDENCIES 
The tool depends on the following: 

    * vmd-python by Robin Betz: https://github.com/Eigenstate/vmd-python
    * Python 3 
    * Numpy
    * Matplotlib and seaborn

## GOALS
Goals here are the following: 
  1. Retrieve PDB file 
  1. Do the edits in the PDB file:
      Remove HETATM HOH molecules 
      Convert ATP to LiG
      Remove ADP and EDD
      Remove CONECT 
  1. Then open up VMD to visualize the given PDB file 
  1. Run the script "sourcce command" in the VMD tk Console. 
  1. Finally, graph the output
  
I'm doing all of this in a single code file which I intend to change in the future, but the 
code does everything in one go. 
