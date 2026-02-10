# mers-md-simulation
Molecular Dynamics Simulation of MERS-CoV virus using OpenMM to find target residue pairs to inhibit the infection. 

The raw simulation results are parsed and analysed to find target residue pairs by,
  1. Calculating the salt bridges between the interchain residues of virus and human receptors.
  2. Calculating the hydrogen bonds occupancies for residues of virus and human receptors.
  3. Calculating the Root-mean-square deviation (RMSD) for all simulations.
  4. Prediction of binding affinity using the energy files of each simulation run.


### 1. Salt Bridge Calculation

Use VMD application to find the salt bridges,
1. Combine all trajectory files in dcd to one dcd file using vmd/dcd_combine.tcl file. (run ```source dcd_combine.tcl``` on vmd terminal)
2. Load the molecule to VMD (raw input pdb file with heteroatoms removed)
3. Load the combined dcd file over the raw file (load data into molecule)
4. Run: Extensions -> Analyis -> Salt Bridges -> Find salt bridges

### 2. Hydrogen bond occupancy rate

On VMD with molecule and combined dcd file loaded as in previous step,

Run: Extensions -> Hydrogen Bonds (check config below) -> select write output to files -> Find hydrogen bonds

Hydrogen Bonds Config (can vary for each):

selection 1: protein and chain A 
selection 2: protein and chain B

Selection 1 is the: Both (donor & acceptor)

Dono-Acceptor Distance (A): 4.0
Angle Cutoff (degrees): 20

calculate detailed info for: Residue_pairs


Note: RMSD calculation and binding affinity prediction steps are part of our another tool under our research lab. The link will be added here as soon as the respective author hosts the code. 
For more information please refer to our paper at: 
