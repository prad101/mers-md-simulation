# MERS-CoV Molecular Dynamics Simulation

This project utilizes **OpenMM** to perform Molecular Dynamics (MD) simulations of the MERS-CoV virus. The primary objective is to identify target residue pairs that could inhibit infection by analyzing the interactions between viral proteins and human receptors.

## Analysis Pipeline
The raw simulation results are processed through four key analytical steps to identify critical residue pairs:
1.  **Salt Bridge Analysis:** Identification of interchain electrostatic interactions.
2.  **Hydrogen Bond Occupancy:** Calculation of bond persistence between the virus and human receptors.
3.  **RMSD Calculation:** Tracking the structural stability of the system over time.
4.  **Binding Affinity Prediction:** Estimating the strength of the interaction using simulation energy files.

---

## 🛠 Analysis Workflow

### 1. Salt Bridge Calculation
Use the **VMD (Visual Molecular Dynamics)** application to identify salt bridges:

1.  **Consolidate Trajectories:** Combine multiple `.dcd` trajectory files into a single file using the `vmd/dcd_combine.tcl` script. Execute `source dcd_combine.tcl` within the VMD terminal.
2.  **Load Molecule:** Load the raw input PDB file, ensuring all heteroatoms have been removed.
3.  **Import Data:** Load the combined `.dcd` file onto the previously loaded molecule.
4.  **Execute Analysis:** Navigate to `Extensions` -> `Analysis` -> `Salt Bridges` -> `Find salt bridges`.

> **Note:** For downstream analysis of salt bridge residual pairs, refer to `sim_analysis_viz.ipynb` notebook.

### 2. Hydrogen Bond Occupancy
With the molecule and consolidated trajectory loaded in VMD:

1.  **Navigate:** Go to `Extensions` -> `Hydrogen Bonds`.
2.  **Configuration:**
    * **Selection 1:** `protein and chain A`
    * **Selection 2:** `protein and chain B`
    * **Selection 1 Role:** Set to `Both (donor & acceptor)`
    * **Distance Cutoff:** 4.0 Å
    * **Angle Cutoff:** 20°
3.  **Output:** Select "Write output to files" and ensure the calculation is set for `Residue_pairs`.

> **Note:** Further visualization and occupancy rate analysis can be found in `sim_analysis_viz.ipynb`.

---

## 📊 Results & Data Structure

All simulation outputs and charts are organized in the `results/` directory:
* **Raw Simulation Data:** `results/simulation_{1/2/3}`
* **Visualizations:** `results/chart`

### Dataset Scale
Each simulation covers a **total duration of 1000 ns**. The data is split into 100 DCD and PDB files per simulation, with each file representing a consecutive 10 ns segment.

---

## 📝 Ongoing Research
The **RMSD calculation** and **binding affinity prediction** tools are currently part of another tool under our research lab. Links to these modules will be added to this repository as soon as they are made public by the respective authors.