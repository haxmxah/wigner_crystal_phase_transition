# Standard library
import os
from pathlib import Path

# Third-party scientific / numeric
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit 

# Visualization
import matplotlib.pyplot as plt
import plotly.express as px
from matplotlib.ticker import LogLocator
# Utils
import warnings

# Export engines
import kaleido
from matplotlib.ticker import MaxNLocator

# from plot_style import set_plot_style

# set_plot_style()

plt.style.use('science.mplstyle')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "STIXGeneral"
})


def get_simulation_parameters(filepath):
    params = {}
    keys = [
        "N", "freeze_mc_steps", "alpha", "density", 
        "charge", "initial_temp", "final_temp", "temp_step"
    ]
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i < len(keys):
                value = line.split('!')[0].strip()
                params[keys[i]] = float(value) if i > 0 else int(value)
                
    return params

def get_simulation_energy_evolution(simulation_path):
    energy_evolution_file = simulation_path / "energy.out"

    df = pd.read_csv(energy_evolution_file, sep='\s+', header=None, names=['mc_step', 'energy', 'temperature', ])
    reset_indices = df['mc_step'].diff() < 0
    offsets = np.where(reset_indices, df['mc_step'].shift(1), 0)
    cumulative_offset = offsets.cumsum()
    df['total_mc_step'] = df['mc_step'] + cumulative_offset

    plt.figure()
    plt.plot(df['total_mc_step'], df['energy'], label='Energía del sistema')
    plt.xlabel('MC step')
    plt.ylabel('Total Energy')
    plt.title('Energy evolution')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(simulation_path / "energy_vs_mc_steps.pdf")
    
    
def get_energy_cv_gamma_evolution_vs_temperature(simulation_path, parameters):
    energy_evolution_file = simulation_path / "heat_capacity.out"
    N = parameters['N']
    density = parameters['density']
    df = pd.read_csv(energy_evolution_file, sep='\s+', header=None, names=['temperature', 'energy', 'cv', 'gamma', 'acc_ratio'])

    # Detection of the phase transition
    max_cv = df['cv'].max()
    Tc = df[df['cv'] == max_cv][['temperature', 'gamma']].values
    max_cv = df['cv'].max()
    Tc_val = Tc[0][0]
    gamma_val = Tc[0][1]
    cv_norm = max_cv / N

    results_dict = {
        'Tc': [Tc_val],
        'gamma': [gamma_val],
        'cv_max_norm': [cv_norm],
        'density' : [density],
        'simulation': [simulation_path],
    }

    df_results = pd.DataFrame(results_dict)

    output_csv = simulation_path / "simulation_results.csv"
    if not os.path.isfile(output_csv):
        df_results.to_csv(output_csv, index=False)
    else:
        df_results.to_csv(output_csv, mode='a', header=False, index=False)

    print(f"Results saved in: {output_csv}")

    # E/N vs. T
    plt.figure()
    plt.plot(df['temperature'], df['energy']/N, label='Energía del sistema')
    plt.xlabel(r'$T$')
    plt.ylabel(r'$E/N$')
    plt.grid(True)
    plt.tight_layout()
    plt.xscale('log')
    plt.savefig(simulation_path / "energy_per_particle_vs_temperature.pdf")
    

    # Cv vs. T
    plt.figure()
    plt.plot(df['temperature'], df['cv']/N, color = 'red', label='Energía del sistema')
    plt.xlabel(r'$T$')
    plt.ylabel(r' $C_v/N$')
    plt.grid(True)
    plt.tight_layout()
    plt.xscale('log')
    plt.savefig(simulation_path / "cv_per_particle_vs_temperature.pdf")
    

    # Gamma vs. T
    plt.figure()
    plt.plot(df['temperature'], df['gamma']/N, color = 'green', label='Energía del sistema')
    plt.xlabel(r'$T$')
    plt.ylabel(r'$\Gamma$')
    plt.grid(True)
    plt.tight_layout()
    plt.xscale('log')
    plt.savefig(simulation_path / "gamma_vs_temperature.pdf")
    

    # Gamma vs. T
    plt.figure()
    plt.plot(df['temperature'], df['acc_ratio'], color = 'purple', label='Energía del sistema')
    plt.xlabel(r'$T$')
    plt.ylabel(r'$\Gamma$')
    plt.grid(True)
    plt.tight_layout()
    plt.xscale('log')
    plt.savefig(simulation_path / "gamma_vs_temperature.pdf")
    

def get_gdr(simulation_path, parameters):
    gdr_file = simulation_path / "rdf.out"
    df = pd.read_csv(gdr_file, sep='\s+', header=None, names=['r', 'gr',])
    rho = parameters['density']
    r = df['r'].to_numpy()
    gr = df['gr'].to_numpy()
    dr = r[1]-r[0]
    
    dNc = 4 * np.pi * (r**2) * rho * gr * dr
    Nc = np.cumsum(dNc)

    peak_idx, _ = find_peaks(gr, height=1.0)

    results = []
    for i, idx in enumerate(peak_idx):
        distance_r = r[idx]
        idx_end_peak = min(idx + 2, len(r) - 1) 
        acc_neigh = Nc[idx_end_peak]
        results.append({
            'r' : distance_r,
            'Acc. Neigh' : int(round(acc_neigh))
        })

    df_results = pd.DataFrame(results)
    output_csv = simulation_path / "gr_peaks_results.csv"

    if not os.path.isfile(output_csv):
        df_results.to_csv(output_csv, index=False)
    else:
        df_results.to_csv(output_csv, mode='a', header=False, index=False)

    print(f"Results saved in: {output_csv}")

    # Plot 
    plt.figure()
    plt.plot(df['r'], df['gr'], color = 'red', label='Energía del sistema')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$g(r)$')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(simulation_path / "gdr.pdf")
    

    return df_results

def get_structure_factor(simulation_path, parameters):
    gdr_file = simulation_path / "sq.out"
    df = pd.read_csv(gdr_file, sep='\s+', header=None, names=['q', 'sq',])

    # Plot 
    plt.figure()
    plt.plot(df['q'], df['sq'], color = 'green', label='Energía del sistema')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$g(r)$')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(simulation_path / "sq.pdf")
    

    return 

def get_positions_plot(filename, saved_filename):

    with open(filename, 'r') as f:
        lines = f.readlines()
    
    num_atoms = int(lines[0])
    block_size = num_atoms + 2
    total_frames = len(lines) // block_size
    
    all_data = []
    for f in range(total_frames):
        start = f * block_size + 2
        end = start + num_atoms
        frame_coords = [l.split()[1:4] for l in lines[start:end]]
        
        df_frame = pd.DataFrame(frame_coords, columns=['x', 'y', 'z'], dtype=float)
        df_frame['frame'] = f  
        all_data.append(df_frame)

    df = pd.concat(all_data, ignore_index=True)

    last_frame = df[df['frame'] == df['frame'].max()]

    fig = px.scatter_3d(
        last_frame,
        x='x', y='y', z='z',
        opacity=1.0,
        template="plotly_white"
    )

    fig.update_traces(
        marker=dict(
            size=7,
            color='red',
            line=dict(width=0.5, color='black')
        )
    )

    # Box edges 
    xmin, ymin, zmin = last_frame[['x','y','z']].min()
    xmax, ymax, zmax = last_frame[['x','y','z']].max()

    edges = [
        ([xmin,xmax],[ymin,ymin],[zmin,zmin]),
        ([xmin,xmax],[ymax,ymax],[zmin,zmin]),
        ([xmin,xmax],[ymin,ymin],[zmax,zmax]),
        ([xmin,xmax],[ymax,ymax],[zmax,zmax]),

        ([xmin,xmin],[ymin,ymax],[zmin,zmin]),
        ([xmax,xmax],[ymin,ymax],[zmin,zmin]),
        ([xmin,xmin],[ymin,ymax],[zmax,zmax]),
        ([xmax,xmax],[ymin,ymax],[zmax,zmax]),

        ([xmin,xmin],[ymin,ymin],[zmin,zmax]),
        ([xmax,xmax],[ymin,ymin],[zmin,zmax]),
        ([xmin,xmin],[ymax,ymax],[zmin,zmax]),
        ([xmax,xmax],[ymax,ymax],[zmin,zmax])
    ]

    for ex, ey, ez in edges:
        fig.add_scatter3d(x=ex, y=ey, z=ez, mode='lines',
                        line=dict(color='black', width=4),
                        showlegend=False)

    fig.update_layout(
        scene=dict(
            aspectmode='cube',
            xaxis=dict(showbackground=False, title='x'),
            yaxis=dict(showbackground=False, title='y'),
            zaxis=dict(showbackground=False, title='z'),
            camera=dict(
                projection=dict(type="orthographic"),
                eye=dict(x=1.6, y=1.6, z=1.3)
            )
        ),
        title_x=0.5
    )

    # Positions array
    pos = last_frame[['x','y','z']].values

    # Bonds
    tree = cKDTree(pos)
    dists, _ = tree.query(pos, k=7)  
    a = np.median(dists[:,1])        
    pairs = tree.query_pairs(r=1.2*a)

    for i, j in pairs:
        fig.add_scatter3d(
            x=[pos[i,0], pos[j,0]],
            y=[pos[i,1], pos[j,1]],
            z=[pos[i,2], pos[j,2]],
            mode='lines',
            line=dict(color='black', width=1),
            opacity=0.5,
            showlegend=False
        )

    fig.update_layout(
        margin=dict(l=0, r=0, b=0, t=0),
        scene=dict(
            aspectmode='cube',
            xaxis=dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, title=''),
            yaxis=dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, title=''),
            zaxis=dict(showbackground=False, showticklabels=False, showgrid=False, zeroline=False, title=''),
            camera=dict(
                projection=dict(type="orthographic"),
                eye=dict(x=1.6, y=1.6, z=1.3)
            ),
            
        )
    )

    fig.write_html(f"{saved_filename}.html")
    fig.write_image(f"{saved_filename}.pdf", scale=2, width=600, height=600)

    return 

from pathlib import Path


def get_phase_transition_diagram(output_dir):
    output_dir = Path(output_dir)

    data = []
    missing = 0

    for folder in sorted(output_dir.iterdir()):
        if folder.is_dir():

            file1 = folder / "simulation_results.csv"
            file2 = folder / "result_simulation.csv"

            if file1.exists():
                file = file1
            elif file2.exists():
                file = file2
            else:
                missing += 1
                warnings.warn(f"No CSV found in: {folder}")
                continue

            df = pd.read_csv(file)

            if {"density", "Tc"}.issubset(df.columns):
                density = df["density"].iloc[0]
                Tc = df["Tc"].iloc[0]
            else:
                raise ValueError(f"Missing columns in {file}")

            data.append((density, Tc))
    
    df = pd.DataFrame(data, columns=["density", "Tc"])
    df = df.sort_values("density")

    print(df)

    # Fit 
    fit_func = lambda Tc, k: k * (Tc**3)
    popt, pcov = curve_fit(fit_func, df["Tc"], df["density"])
    k = popt[0]
    Tc_vals = np.linspace(df["Tc"].min(), df["Tc"].max(), 100)
    rho_theoretical = fit_func(Tc_vals, k)

    # Plot 
    plt.figure()

    plt.fill_between(Tc_vals, df["density"].min() -0.1, rho_theoretical, 
                    color='lightgreen', alpha=0.3)
    plt.fill_between(Tc_vals, rho_theoretical, df["density"].max(), 
                    color='salmon', alpha=0.3)

    plt.plot(Tc_vals, rho_theoretical, "--", color='black', lw=1, 
            label=rf"$\rho = {k:.2f} T_c^3$")
    plt.plot(df["Tc"], df["density"], "o", markersize=4, 
            label="MC Simulation")

    plt.text(df["Tc"].min() + 0.9, df["density"].min()+0.001, "Gas", 
            fontsize=12, color='darkgreen')
    plt.text(df["Tc"].max() - 0.5, df["density"].max() - 800, "Coulomb Crystal", 
            fontsize=12, color='darkred', horizontalalignment='right')

    ymin = df["density"].min()
    ymax = df["density"].max()

    plt.ylim(ymin, ymax)

    plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=6))
    plt.xlabel("$T_c$")
    plt.ylabel(r"$\rho$")
    plt.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.20),
           ncol=2) 

    # plt.grid(True)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig( output_dir / 'phase_transition_diagram.pdf', bbox_inches='tight')
    
    print(f"Phase diagram generated. K = {k:.2f}")

    
    return

if '__main__':
    
    output_dir = Path("output_n216")

    # for simulation_path in sorted(output_dir.iterdir()):
    #     if not simulation_path.is_dir():
    #         continue

    #     try:
    #         print(f"\nProcessing: {simulation_path}")

    #         input_parameters_file = simulation_path / "input_parameters.in"
    #         filename = simulation_path / "final_position.xyz"
    #         saved_filename = simulation_path / "positions"

    #         if not input_parameters_file.exists():
    #             warnings.warn(f"Missing input_parameters.in in {simulation_path}")
    #             continue

    #         parameters = get_simulation_parameters(input_parameters_file)

    #         get_simulation_energy_evolution(simulation_path)
    #         get_energy_cv_gamma_evolution_vs_temperature(simulation_path, parameters)
    #         get_positions_plot(filename, saved_filename)
    #         get_gdr(simulation_path, parameters)

    #     except Exception as e:
    #         warnings.warn(f"Error processing {simulation_path}: {e}")

    get_phase_transition_diagram(output_dir)
