
from openmm.app import PDBFile, Modeller
import pdbfixer
import math
import os
import pandas as pd
import numpy as np
import seaborn as sns
from plotnine import *
from matplotlib import pyplot as plt
import math
from matplotlib.lines import Line2D # Import Line2D for the legend

def cif_to_pdb(name):
    # cif to pdb without pdbfixer
    path = f'input/{name}.cif'
    cif = PDBFile(path)
    topology = cif.getTopology()
    positions = cif.getPositions()

    # Write the fixed structure to a temporary PDB file
    with open(f'input/{name}.pdb', 'w') as f:
        PDBFile.writeFile(topology, positions, f, keepIds=True)


def cif_to_pdbfix(name):
    # cif to pdb with pdbfixer
    path = f'input/{name}.cif'
    fixer = pdbfixer.PDBFixer(filename=path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    # Write the fixed structure to a temporary PDB file
    with open(f'input/{name}.pdb', 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)


def filter_atoms_and_save(pdb):    
    topology = pdb.topology
    positions = pdb.positions
    modeller = Modeller(topology, positions)

    hetero_atoms = [r for r in topology.residues() if (r.name == 'NAG') or (r.name == 'BMA')]
    modeller.delete(hetero_atoms)

    # with open('input/filter_4l.pdb', 'w') as f:
    #     PDBFile.writeFile(modeller.topology, modeller.positions, f)

    return modeller


def inter_chain_selection(sb_path, sb_path_new):
    os.makedirs(sb_path_new, exist_ok=True)
    for file in os.listdir(sb_path):
        chain1 = file.split('_')[1].split('-')[0]
        chain2 = file.split('_')[2].split('.')[0]

        if chain1 != chain2:
            with open(f'{sb_path}/{file}', 'r') as f:
                content = f.read()

            with open(f'{sb_path_new}/{file}', 'w') as f:
                f.write(content)
                print(f"Inter chain salt bridges save in, {sb_path_new}")


def sb_extract_residues(residue_pair):
   parts = residue_pair.split('-')
   
   if 'chainA' in parts[0]:
       chainA_residue = parts[0].split('_')[0]  
       chainB_residue = parts[1].split('_')[0]  
   else:  
       chainA_residue = parts[1].split('_')[0]  
       chainB_residue = parts[0].split('_')[0]  
   
   return pd.Series([chainA_residue, chainB_residue])

# salt bridge plot (cummulative)
def plot_salt_bridge_grid_charts(df: pd.DataFrame, simulation_time: int):

    save_dir = "charts/"
    unique_bridges = df['sb_name'].unique()
    n_bridges = len(unique_bridges)
    
    rows = 4
    cols = 2
    rows = math.ceil(n_bridges / cols)

    panel_labels = 'ABCDEFG'[:n_bridges]
    unique_funnels = df['funnel'].unique()
    colors = ['tab:blue', 'tab:orange', 'tab:green']

    # Setup the plot
    sns.set(style="whitegrid")
    palette = dict(zip(unique_funnels, colors))

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True, sharey=True)
    axes = axes.flatten()
    plt.setp(axes, xlabel='Time (seconds)')

    # Plotting each salt bridge in a separate subplot
    for i, sb in enumerate(unique_bridges):
        ax = axes[i]
        sub_df = df[df['sb_name'] == sb]
        virus_human_pair = sub_df['virus_human'].iloc[0]

        data_thresh = 0.25 * max(sub_df['timestamp'])
        sub_df = sub_df[sub_df['timestamp'] >= data_thresh]
        # ymin = math.floor(sub_df['sb_val'].min())
        # ymax = math.ceil(sub_df['sb_val'].max())

        sns.lineplot(
            data=sub_df,
            x="timestamp",
            y="sb_val",
            hue="funnel",
            ax=ax,
            palette=palette,
            linewidth=1.5,
            legend=False 
            
        )
        ax.axhline(4, color='black', linestyle='--', linewidth=0.8)
        ax.set_title(f"{virus_human_pair}", fontweight='bold', fontsize=14)
        ax.set_ylabel("Distance (Å)", fontweight='bold', fontsize=15)
        ax.set_xlabel("Time (ns)", fontweight='bold', fontsize=14, )

        ax.text(0.0, 1.1, panel_labels[i], transform=ax.transAxes, 
                fontsize=16, fontweight='bold', va='top', ha='right', )
        # ax.set_ylim(ymin, ymax)
        # ax.legend().set_title("")
        ax.tick_params(labelbottom=True)
        last_ax = ax

    if last_ax:
        legend_ax = axes[n_bridges] # Use the next available slot
        legend_ax.axis('off') # Turn off the plot box

        # Get handles and labels from the last plot
        handles, labels = last_ax.get_legend_handles_labels()
        
        # Manually add the dashed line to the legend
        blue_line = Line2D([0], [0], color='tab:blue', linestyle='-',
                            linewidth=3.0, label='Simulation 1')
        orange_line = Line2D([0], [0], color='tab:orange', linestyle='-', 
                            linewidth=3.0, label='Simulation 2')
        green_line = Line2D([0], [0], color='tab:green', linestyle='-', 
                            linewidth=3.0, label='Simulation 3')
        dashed_line = Line2D([0], [0], color='black', linestyle='--', 
                            linewidth=3.0, label='Interaction Cutoff (4.0 Å)')
    
        handles.append(blue_line)
        handles.append(orange_line)
        handles.append(green_line)
        handles.append(dashed_line)

        labels.append(f'Simulation 1 ({simulation_time} ns)')
        labels.append(f'Simulation 2 ({simulation_time} ns)')
        labels.append(f'Simulation 3 ({simulation_time} ns)')
        labels.append('Interaction Cutoff (4.0 Å)')

        # Create the legend in the new axis
        legend_ax.legend(handles, labels, loc='center', fontsize=14, 
                        title_fontsize=16, frameon=True, shadow=True, title='Legend')

    # Remove any unused axes
    for j in range(n_bridges + 1, len(axes)):
        fig.delaxes(axes[j])
    

    plt.tight_layout()
    plt.suptitle("Salt Bridge Analysis between Mers-CoV and DPP-4", fontsize=18, fontweight='bold',  y=1.02, fontfamily='sans-serif')
    plt.savefig(f"{save_dir}/salt_bridge_all.png", dpi=350)
    plt.show()

# converts frame to timestamp
def frame_to_timestamp(df, step_size):

    # converts to ns - 2 * (step_size / 1,000,000)
    df['timestamp'] = df['frames'] * (2 * step_size)/1_000_000
    return df


def salt_bridge_df(sb_path_new, step_size, funnel_name):
    df_list = []
    for file in os.listdir(sb_path_new):
        df = pd.read_csv(f'{sb_path_new}/{file}', delimiter=' ', header=None)
        df.rename({0:'frames', 1: 'sb_val'}, axis=1, inplace=True)
        sb_name = str(file).split(".")[0]
        df['sb_name'] = sb_name.split("saltbr-")[1]
        df = frame_to_timestamp(df, step_size)
        df['funnel'] = funnel_name
        df_list.append(df)

    sb_df = pd.concat(df_list)
    return sb_df


def h2_load_viz(h2_path, step_size):
    
    df = pd.read_csv(h2_path, delimiter=' ', header=None)
    df.rename({0:'frames', 1: 'h_bonds'}, axis=1, inplace=True)
    df = frame_to_timestamp(df, step_size)

    plot = (
        ggplot(df)
        + aes(x="timestamp", y="h_bonds")
        + labs(
            x="Simulation time (ns)",
            y="No. of H-bonds",
            title="H-bonds comparison",
        )
        + geom_line()
        + theme(figure_size=(8, 6), 
                title = element_text(size = 8),
                text = element_text(size = 6))
    # + facet_grid("sb_name ~ .", scales="free_y")
    )
    plot.show()

# h2 file preprocessor
def h2_detail_load_preprocess(h2_detail_path, simulation_num):
    
    df = pd.read_csv(h2_detail_path, delimiter='\t ', header=1)
    
    df['occupancy'] = df['occupancy'].apply(lambda x: float(str(x).strip('%')))
    df.rename({'donor \t\t': 'donor', 'acceptor \t': 'acceptor'}, axis=1, inplace=True) 

    df['donor_chain_B'] = df.apply(lambda x: x['donor'].split("-")[1] 
                                                       if x['donor'].startswith('SegB') else x['acceptor'].split("-")[1], 
                                                       axis=1)
    df['acceptor_chain_A'] = df.apply(lambda x: x['acceptor'].split("-")[1] 
                                                          if x['acceptor'].startswith('SegA') else x['donor'].split("-")[1], 
                                                          axis=1)
    
    # renumbering the amino acids to original numbering
    df['acceptor_chain_A'] = df['acceptor_chain_A'].apply(lambda x: x[:3] + (str(39 + (int(x[3:])) - 1)))
    df['donor_chain_B'] = df['donor_chain_B'].apply(lambda x: x[:3] + (str(382 + (int(x[3:])) - 1)))

    df['Donor - Acceptor'] = df['donor_chain_B'] + ":" + df['acceptor_chain_A']

    df['funnel'] = f'Simulation {simulation_num}'

    return df

# h2 occupancy plot (individual)
def h2_occupancy_viz(h2_detail_df, sim_num):

    save_dir = "charts/h2"
    os.makedirs(save_dir, exist_ok=True)
    occupancy_theshold = 60

    h2_detail_df = h2_detail_df.loc[h2_detail_df['occupancy'] > occupancy_theshold].reset_index(drop=True)
    h2_detail_df = h2_detail_df.sort_values('occupancy', ascending=False).reset_index(drop=True)

    color_palette = ['white', 'tab:blue', 'tab:orange', 'tab:green']
    plt.figure(figsize=(12, 8))
    order = h2_detail_df['Donor - Acceptor']

    ax = sns.barplot(data=h2_detail_df, y='Donor - Acceptor', x='occupancy', order=order, color=color_palette[sim_num], errorbar=None)
    ax.bar_label(ax.containers[0], fontsize=14, fontweight='bold', fontfamily='sans-serif')
    ax.set_xlim(left=50)

    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(14)
        label.set_fontfamily('sans-serif')

    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(14)
        label.set_fontfamily('sans-serif')

    plt.xlabel('Occupancy (%)', fontsize=17, fontweight='bold', fontfamily='sans-serif')
    plt.ylabel('Virus - Human', fontsize=17, fontweight='bold', fontfamily='sans-serif')
    plt.title(f'Mers-CoV and DPP-4 Hydrogen bonds - Simulation {sim_num} \n(>{occupancy_theshold}% occupancy, {h2_detail_df.shape[0]} pairs)', fontsize=17, fontweight='bold', fontfamily='sans-serif')
    plt.tight_layout()
    plt.savefig(f"{save_dir}/simulation_{sim_num}.png", dpi=300)
    plt.show()

# h2 occupancy average over simulation plot (cummulative)
def h2_occupancy_viz_overlap_avg(h2_detail_df):

    save_dir = "charts/h2"
    os.makedirs(save_dir, exist_ok=True)

    color_palette = ['tab:blue', 'tab:orange', 'tab:green']
    plt.figure(figsize=(12, 8))

    h2_detail_df["Donor - Acceptor"] = h2_detail_df["Donor - Acceptor"].apply(lambda x : x.replace(":", ":\n"))
    h2_detail_df = h2_detail_df.sort_values(by=['funnel', 'occupancy'], ascending=[True, False]).reset_index(drop=True)
    # labels = [x.replace(":", ":\n") for x in labels] 

    top_5 = h2_detail_df.groupby("Donor - Acceptor")["occupancy"].mean().nlargest(5).index

    graph_data = h2_detail_df[h2_detail_df["Donor - Acceptor"].isin(top_5)]
    # order = graph_data.sort_values('occupancy', ascending=False)['Donor - Acceptor']

    ax = sns.barplot(data=graph_data
            , y='Donor - Acceptor', x='occupancy', hue='funnel', order=top_5, errorbar=None)

    for container in ax.containers:
        ax.bar_label(container, fontsize=14, fontweight='bold', fontfamily='sans-serif')

    ax.set_xlim(0, 100)

    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(14)
        label.set_fontfamily('sans-serif')

    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(14)
        label.set_fontfamily('sans-serif')

    legend = ax.get_legend()
    legend.set_bbox_to_anchor((1.05, 1.0))
    # legend._set_loc('upper left')
    legend.set_title("Funnel", prop={"size": 14, "family": "sans-serif", "weight": "bold"})

    # legend._set_loc('upper left')

    if legend is not None:
        for text in legend.get_texts():
            # text.set_fontweight('bold')
            text.set_fontsize(14)
            text.set_fontfamily('sans-serif')

    # for tick in ax.get_yticks():
    #     tick.set_fontweight('bold')

    plt.xlabel('Occupancy (%)', fontsize=17, fontweight='bold', fontfamily='sans-serif')
    plt.ylabel('Virus - Human', fontsize=17, fontweight='bold', fontfamily='sans-serif')
    plt.title(f'Mers-CoV and DPP-4 Hydrogen bonds - Overlaps across simulations \n(Top {len(top_5)} pairs by Average Occupancy)', fontsize=17, fontweight='bold', fontfamily='sans-serif')
    plt.tight_layout()
    plt.savefig(f"{save_dir}/simulation_average.png", dpi=300)
    plt.show()
