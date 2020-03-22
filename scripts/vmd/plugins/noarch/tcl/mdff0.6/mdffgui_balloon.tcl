#text for the popup balloons when user hover overs GUI elements

proc MDFFGUI::gui::workingdir {} {
  return "Directory to where all output files will be written."
}

proc MDFFGUI::gui::currentmol {} {
  return "Already loaded molecule with a .psf and .pdb for your starting structure."
}

proc MDFFGUI::gui::currentdensitymap {} {
  return "Add density map(s) as a target for fitting."
}

proc MDFFGUI::gui::fixedpdbtext {} {
  return "atomselection text for fixing a part of your structure during fitting."
}

proc MDFFGUI::gui::fixedpdbcolumn {} {
  return "Column in output PDB file to use to mark fixed atoms."
}

proc MDFFGUI::gui::simoutputname {} {
  return "Name prefix for simulation output files."
}

proc MDFFGUI::gui::simstep {} {
  return "Step of the fitting process. Numbers after 1 will automatically proceed from where the previous simulation finished."
}

proc MDFFGUI::gui::restraints {} {
  return "Automatically generate optional simulation restraints to prevent overfitting."
}

proc MDFFGUI::gui::ssrestraintstip {} {
  return "Restrain dihedral angles in secondary structures of protein as well as dihedral angles and base pairs in nucleic acids."
}

proc MDFFGUI::gui::chiralityrestraints {} {
  return "Restrain chiral centers to their initial handedness."
}

proc MDFFGUI::gui::cispeptiderestraints {} {
  return "Restrain cis peptide bonds to their initial cis or trans configuration."
}

proc MDFFGUI::gui::parameters {} {
  return "Force field parameter files to use during simulation."
}

proc MDFFGUI::gui::mdfffiles {} {
  return "Input and output file information for MDFF fitting simulations."
}

proc MDFFGUI::gui::simparams {} {
  return "Parameters controlling length, environment, and other general aspects of the simulation."
}

proc MDFFGUI::gui::itemp {} {
  return "Initial temperature of the simulation in Kelvin."
}

proc MDFFGUI::gui::ftemp {} {
  return "Final temperature of the simulation in Kelvin."
}

proc MDFFGUI::gui::minsteps {} {
  return "Number of minimization steps to perform prior to simulation."
}

proc MDFFGUI::gui::numsteps {} {
  return "Number of simulation steps to perform."
}

proc MDFFGUI::gui::gridforcelite {} {
  return "Perform simpler trilinear interpolation of fitting forces. Faster performance with slightly less accuracy."
}

proc MDFFGUI::gui::gridforceoff {} {
  return "Turn fitting forces off."
}

proc MDFFGUI::gui::pbc {} {
  return "Use periodic boundary conditions (for use with explicit solvent environments)."
}

proc MDFFGUI::gui::gbis {} {
  return "Use Generalized Born Implicit Solvent."
}

proc MDFFGUI::gui::vac {} {
  return "Perform simulation in vacuum."
}

proc MDFFGUI::gui::xmdff {} {
  return "Parameters for xMDFF (MDFF for low-resolution X-ray crystallography)."
}

proc MDFFGUI::gui::doxmdff {} {
  return "Setup xMDFF simulation (MDFF for low-resolution X-ray crystallography)."
}

proc MDFFGUI::gui::xmdffdensity {} {
  return "Reflection data and parameters for generating real-space density map from X-ray diffraction data."
}

proc MDFFGUI::gui::xmdffsel {} {
  return "atomselection text for part of structure to refine."
}

proc MDFFGUI::gui::xmdffsteps {} {
  return "Number of refinement steps in each iteration of xMDFF (i.e., how many fitting steps to perform before re-generating real-space density map."
}

proc MDFFGUI::gui::xmdffsymmetry {} {
  return "PDB file with properly formatted CRYST line for setting symmetry group of diffraction data. Only needed if symmetry group is not already in input diffraction data."
}

proc MDFFGUI::gui::xmdffbf {} {
  return "Calculate beta factors during each map generation step."
}

proc MDFFGUI::gui::imd {} {
  return "Interactive MDFF to manually manipulate the structure during fitting."
}

proc MDFFGUI::gui::imdport {} {
  return "Port to use for IMD to communicate between VMD and NAMD."
}

proc MDFFGUI::gui::imdfreq {} {
  return "Frequency of updating VMD with NAMD timesteps."
}

proc MDFFGUI::gui::imdkeep {} {
  return "Save every X number of NAMD timesteps as frames in VMD."
}

proc MDFFGUI::gui::doimd {} {
  return "Turn on interactive MDFF to manually manipulate the structure during fitting."
}

proc MDFFGUI::gui::imdignore {} {
  return "Ignore IMD forces during simulation. Use this option to monitor and analyze fitting progress without manually manipulating the structure."
}

proc MDFFGUI::gui::imdhost {} {
  return "Server running NAMD."
}

proc MDFFGUI::gui::imdwait {} {
  return "Do not proceed with NAMD simulation until VMD has connected with NAMD after startup."
}

proc MDFFGUI::gui::imdproc {} {
  return "Number of CPU cores to run NAMD with."
}
