from openmm.app import *
from openmm import *
from openmm.unit import *
import argparse
from sys import stdout
import os
from tqdm import tqdm
from utils import cif_to_pdb, cif_to_pdbfix, filter_atoms_and_save
from datetime import datetime

# Set the CUDA platform to enable GPU acceleration
platform = Platform.getPlatformByName('CUDA')

# Specify the device index to use all available GPUs
properties = {'DeviceIndex': '0,1,2,3'} # Assuming 4 GPUs are indexed from 0 to 3

def run_simulations(name, n_femto, simulation_number):

    # Load the CIF file and fix it using PDBFixer (only if required)
    cif_to_pdbfix(name)
    pdb = PDBFile(f'input/{name}.pdb')
    print(pdb.topology, n_femto)

    # Create a system using a force field
    forcefield = ForceField('amber14-all.xml', 'tip3pfb.xml') # Adjust the force field files as necessary

    # retain only atoms and add hydrogen
    pdb = filter_atoms_and_save(pdb)

    # Create the system from the fixed PDB file
    system = forcefield.createSystem(pdb.topology, 
                                    nonbondedMethod=PME,
                                    nonbondedCutoff=1*nanometer, 
                                    constraints=HBonds,
                                    rigidWater=True,            
                                    removeCMMotion=True)

    integrator = LangevinIntegrator(310.45*kelvin, 1/picosecond, n_femto*femtoseconds)
    simulation = Simulation(pdb.topology, 
                            system, 
                            integrator, 
                            platform, 
                            properties
                            )

    simulation.context.setPositions(pdb.positions)

    # Output directory for PDB files
    output_pdb_directory = os.path.join(f"results_{name}_{simulation_number}", f"output_pdb_files", 'MT_Tau_Reset_100ns')
    output_dcd_directory = os.path.join(f"results_{name}_{simulation_number}", f"output_dcd_files", 'MT_Tau_Reset_100ns')
    output_energy_directory = os.path.join(f"results_{name}_{simulation_number}", f"output_energy_files", 'MT_Tau_Reset_100ns')
    checkpoint_directory = os.path.join(f"results_{name}_{simulation_number}", f"output_pdb_files", 'MT_Tau_Reset_100ns_checkpoint')
    os.makedirs(output_pdb_directory, exist_ok=True)
    os.makedirs(output_dcd_directory, exist_ok=True)
    os.makedirs(output_energy_directory, exist_ok=True)
    os.makedirs(checkpoint_directory, exist_ok=True)

    # Set up simulation parameters
    # steps = 500500000
    # report_interval = 500000
    # checkpoint_interval = 500000

    # steps = 100100000
    # report_interval = 100000
    # checkpoint_interval = 100000

    # steps = 1000000
    # report_interval = 10000
    # checkpoint_interval = 10000

    # steps = 100000000
    # report_interval = 1000000
    # checkpoint_interval = 1000000

    #for 1000ns simulation time (500M * 2fs = 1us)
    steps = 500_000_000
    report_interval = 5_000_000
    checkpoint_interval = 5_000_000

    # Minimize energy and set initial velocities
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(310.45*kelvin)

    # Check if there are checkpoint files to resume from
    checkpoint_files = [filename for filename in os.listdir(checkpoint_directory) if filename.startswith('checkpoint_')]
    if checkpoint_files:
        # Extract the step number from the checkpoint file names
        step_numbers = [int(''.join(filter(str.isdigit, filename))) for filename in checkpoint_files]
        current_step = max(step_numbers)
        checkpoint_filename = os.path.join(checkpoint_directory, f'checkpoint_{current_step}.chk')
        print(f"Resuming simulation from step {current_step} using {checkpoint_filename}")
        simulation.loadCheckpoint(checkpoint_filename)  # Load the latest checkpoint
    else:
        current_step = 0
        print("Starting a new simulation")

    # Report state data every report_interval steps
    simulation.reporters.append(StateDataReporter(stdout, report_interval, step=True, potentialEnergy=True, temperature=True, volume=True, density=True, progress=True, remainingTime=True, speed=True, totalSteps=steps))

    energy_file_path = os.path.join(output_energy_directory, 'energy.txt')
    energy_file = open(energy_file_path, 'a')

    #dcd_reporter = DCDReporter(os.path.join(output_dcd_directory, 'trajectory.dcd'), report_interval)
    #simulation.reporters.append(dcd_reporter)

    # Save coordinates at specific steps
    for step in tqdm(range(current_step, steps, report_interval)):
        pdb_reporter = PDBReporter(os.path.join(output_pdb_directory, f'output_{step}.pdb'), report_interval)
        pdb_reporter.report(simulation, simulation.context.getState(getPositions=True))

        # DCD reporter
        dcd_reporter = DCDReporter(os.path.join(output_dcd_directory, f'output_{step}.dcd'), report_interval)
        dcd_reporter.report(simulation, simulation.context.getState(getPositions=True))
        
        # Energy reporter
        #energy_file = open(os.path.join(output_energy_directory, f'energy_{step}.txt'), 'w')
        #energy_reporter = StateDataReporter(energy_file, 1, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True)
        
        energy_reporter = StateDataReporter(energy_file, 1, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, progress=True, remainingTime=True, speed=True, totalSteps=steps)

        energy_reporter.report(simulation, simulation.context.getState(getEnergy=True))
        #energy_file.close()
        
        simulation.step(report_interval)

        checkpoint_filename = os.path.join(checkpoint_directory, f'checkpoint_{step+report_interval}.chk')  # Save checkpoint in the subfolder
        simulation.saveCheckpoint(checkpoint_filename)
        print(f"Saved checkpoint at step {step+report_interval} in {checkpoint_directory}")


    energy_file.close()

    # Save final checkpoint
    simulation.saveCheckpoint('final_checkpoint.chk')
    print(f"Simulation completed. Final checkpoint saved.")


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_name", help="Name of protein (if file in input/ is 4l72.pdb, just enter 4l72)")
    parser.add_argument("n_femto", help="Femtoseconds to modify step size (start with 2/3/4)" )
    parser.add_argument("simulation_number", help="simulation_number to differentiate the replication count" )

    args = parser.parse_args()
    pdb_name = str(args.pdb_name)
    n_femto = int(args.n_femto)
    simulation_number = int(args.simulation_number)

    run_simulations(pdb_name, n_femto, simulation_number)
