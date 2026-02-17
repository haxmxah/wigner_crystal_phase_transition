###
# Computes the structure factor S(q) from the simulation output and saves it to a CSV file for later visualization.
# Runs a forran kernel for efficient computation of S(q) across multiple frames of the simulation.
###


import sq_fast
import numpy as np
import pandas as pd
from pathlib import Path

N_MAX = 20
MAX_FRAMES = 5_000

def structure_factor(L, simulation_path):
    n_range = np.arange(-N_MAX, N_MAX + 1)
    nx, ny, nz = np.meshgrid(n_range, n_range, n_range)
    n_vecs = np.vstack([nx.ravel(), ny.ravel(), nz.ravel()]).T
    n_vecs = n_vecs[np.any(n_vecs != 0, axis=1)] # Delete (0,0,0)
    
    q_vecs = (2.0 * np.pi / L) * n_vecs
    q_mags = np.linalg.norm(q_vecs, axis=1)
    
    # Fortran array preparation
    qx = np.ascontiguousarray(q_vecs[:, 0], dtype=np.float64)
    qy = np.ascontiguousarray(q_vecs[:, 1], dtype=np.float64)
    qz = np.ascontiguousarray(q_vecs[:, 2], dtype=np.float64)
    
    q, indices, counts = np.unique(np.round(q_mags, 4), return_inverse=True, return_counts=True)
    sq = np.zeros_like(q)
    
    num_frames = 0

    with open(simulation_path / "positions_t0.xyz", 'r') as f:
        while num_frames < MAX_FRAMES:
            line = f.readline()
            if not line: break
            try:
                N_atoms = int(line.strip())
                f.readline() 
                
                coords_raw = [f.readline().split() for _ in range(N_atoms)]
                coords = np.array(coords_raw)[:, 1:4].astype(np.float64)
                
                x = np.ascontiguousarray(coords[:, 0])
                y = np.ascontiguousarray(coords[:, 1])
                z = np.ascontiguousarray(coords[:, 2])

                sq_frame = sq_fast.sq_calc.compute_sq_kernel(x, y, z, qx, qy, qz)
                
                # Accumulate isotropic S(q)
                s_q_iso_frame = np.zeros_like(q)
                np.add.at(s_q_iso_frame, indices, sq_frame)
                sq += (s_q_iso_frame / counts)
                
                num_frames += 1

            except ValueError:
                break
    
    sq /= num_frames

    df = pd.DataFrame({'q': q, 'S(q)': sq})
    df.to_csv(simulation_path / "sq.csv", index=False)

    return q, sq

if __name__ == "__main__":
    N = 216
    rho = 1.00
    L = (N/ rho)**(1/3) # Can be also read from parameters file if needed
    simulation_path = Path(f"output/n216_density1.00_t0.0001-10-0.1_melted")
    q, sq = structure_factor(L = L, simulation_path=simulation_path)