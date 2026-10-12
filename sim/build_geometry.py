"""
Build and tag the 3D computational domain for LEAF GBC2 pot.
Tags voxels as:
  0: SOLID_WALL (acrylic pot body and lid)
  1: AIR_HEADSPACE (open chamber / plant shoot growth zone)
  2: POROUS_SUBSTRATE (root substrate bed, Y <= 8 mm)
  3: INLET_PORT (Port 1 at Z = -21.6 mm)
  4: OUTLET_PORT (Port 2 at Z = 0.0 mm)
  5: CLOSED_PORT (Port 3 RNAlater at Z = +21.6 mm)
  6: PLANT_CANOPY (leaf surface zone in compartments 1, 2, 3)
"""
import os
import trimesh
import numpy as np
from scipy import ndimage

def build_grid(pitch=0.001): # 1.0 mm pitch by default (fast & robust)
    print(f"Loading CAD assembly with pitch = {pitch*1000:.2f} mm...")
    lid = trimesh.load('cad/lid.step').geometry['M09025F-MDL-2']
    pot = trimesh.load('cad/pot.step').geometry['M09024F-MDL-5']

    # Transformation matrix: rotate 90 deg around Y, translate to mate at pot rim
    T = np.array([
        [ 0.0, 0.0, 1.0, 0.0],
        [ 0.0, 1.0, 0.0, (33.670 - 1.270)/1000.0],
        [-1.0, 0.0, 0.0, 0.0],
        [ 0.0, 0.0, 0.0, 1.0]
    ])
    lid_mated = lid.copy()
    lid_mated.apply_transform(T)
    assembly = trimesh.util.concatenate([pot, lid_mated])

    # Voxelize
    vox = assembly.voxelized(pitch=pitch)
    matrix = vox.matrix.copy() # True = solid wall

    # Coordinate mapping: voxel index (i, j, k) -> world (x, y, z)
    inv_T = np.linalg.inv(vox.transform)
    def world_to_vox(x, y, z):
        p = np.array([x, y, z, 1.0])
        v = inv_T @ p
        return int(round(v[0])), int(round(v[1])), int(round(v[2]))

    def vox_to_world(i, j, k):
        p = np.array([i, j, k, 1.0])
        w = vox.transform @ p
        return w[0], w[1], w[2]

    # Port coordinates
    # Center X = -15.494 mm, Y_top = 38.1 mm
    # Port 1 (Inlet): Z = -21.59 mm
    # Port 2 (Outlet): Z = 0.00 mm
    # Port 3 (RNAlater): Z = +21.59 mm
    r_port = 0.0015875 # 1.5875 mm (dia = 3.175 mm)
    
    ports = {
        'inlet': (-0.015494, -0.02159),
        'outlet': (-0.015494, 0.0),
        'rnalater': (-0.015494, 0.02159)
    }

    # Temporarily seal the top ports to isolate interior cavity from ambient
    capped_matrix = matrix.copy()
    for name, (px, pz) in ports.items():
        vx, vy, vz = world_to_vox(px, 0.036, pz)
        rv = int(np.ceil(r_port / pitch)) + 1
        for di in range(-rv, rv+1):
            for dk in range(-rv, rv+1):
                if di*di + dk*dk <= (rv+1)*(rv+1):
                    ii = np.clip(vx + di, 0, matrix.shape[0]-1)
                    kk = np.clip(vz + dk, 0, matrix.shape[2]-1)
                    capped_matrix[ii, vy-2:min(vy+4, matrix.shape[1]), kk] = True

    # Flood fill exterior ambient space
    padded = np.pad(capped_matrix, 1, mode='constant', constant_values=False)
    labeled, _ = ndimage.label(~padded)
    exterior_mask = (labeled == labeled[0, 0, 0])[1:-1, 1:-1, 1:-1]

    # Interior cavity = not solid and not exterior
    interior = (~matrix) & (~exterior_mask)

    # Classification grid: 0 = solid, 1 = air, 2 = porous
    grid = np.zeros(matrix.shape, dtype=np.int32)
    
    # Classify voxels based on height Y and spatial coordinates
    y_coords = np.array([vox_to_world(0, j, 0)[1] for j in range(matrix.shape[1])])
    z_coords = np.array([vox_to_world(0, 0, k)[2] for k in range(matrix.shape[2])])
    x_coords = np.array([vox_to_world(i, 0, 0)[0] for i in range(matrix.shape[0])])

    # Porous substrate relocated up in Y to chamber bottom shelf: Y = 7.62 mm to 14.62 mm
    j_sub_min = np.searchsorted(y_coords, 0.00762)
    j_sub_max = np.searchsorted(y_coords, 0.01462)
    
    # Split into 3 discrete parts fitting inside the bottom of each chamber along Z:
    # Chamber 1: Z in [-30.5 mm, -9.5 mm]
    # Chamber 2: Z in [-7.5 mm, +7.5 mm]
    # Chamber 3: Z in [+9.5 mm, +30.5 mm]
    # In X: within the plant growth cavities [-11.5 mm, +21.5 mm] (excluding side manifold)
    k_ch1 = (z_coords >= -0.0305) & (z_coords <= -0.0095)
    k_ch2 = (z_coords >= -0.0075) & (z_coords <= +0.0075)
    k_ch3 = (z_coords >= +0.0095) & (z_coords <= +0.0305)
    i_cham = (x_coords >= -0.0115) & (x_coords <= +0.0215)
    
    y_sub_mask = np.zeros(matrix.shape[1], dtype=bool)
    y_sub_mask[j_sub_min:j_sub_max] = True
    
    sub1 = interior & i_cham[:, None, None] & y_sub_mask[None, :, None] & k_ch1[None, None, :]
    sub2 = interior & i_cham[:, None, None] & y_sub_mask[None, :, None] & k_ch2[None, None, :]
    sub3 = interior & i_cham[:, None, None] & y_sub_mask[None, :, None] & k_ch3[None, None, :]
    porous_mask = sub1 | sub2 | sub3

    # Fill interior
    grid[interior] = 1 # Air headspace
    grid[porous_mask] = 2 # Porous substrate (3 partitioned blocks)

    # Tag Ports
    # Tag Inlet (Port 1: 3) and Outlet (Port 2: 4)
    rv_port = int(np.ceil(r_port / pitch))
    
    # Inlet port voxels
    vx1, vy1, vz1 = world_to_vox(ports['inlet'][0], 0.035, ports['inlet'][1])
    for di in range(-rv_port, rv_port+1):
        for dk in range(-rv_port, rv_port+1):
            if di*di + dk*dk <= rv_port*rv_port:
                ii = np.clip(vx1 + di, 0, grid.shape[0]-1)
                kk = np.clip(vz1 + dk, 0, grid.shape[2]-1)
                # Mark upper channel
                for jj in range(max(0, vy1-3), min(grid.shape[1], vy1+4)):
                    if not matrix[ii, jj, kk]:
                        grid[ii, jj, kk] = 3 # INLET

    # Outlet port voxels
    vx2, vy2, vz2 = world_to_vox(ports['outlet'][0], 0.035, ports['outlet'][1])
    for di in range(-rv_port, rv_port+1):
        for dk in range(-rv_port, rv_port+1):
            if di*di + dk*dk <= rv_port*rv_port:
                ii = np.clip(vx2 + di, 0, grid.shape[0]-1)
                kk = np.clip(vz2 + dk, 0, grid.shape[2]-1)
                for jj in range(max(0, vy2-3), min(grid.shape[1], vy2+4)):
                    if not matrix[ii, jj, kk]:
                        grid[ii, jj, kk] = 4 # OUTLET

    # Tag plant canopy leaves seated above the substrate:
    # Chamber 1 & 3: 32 Brassica rapa seeds/plants per chamber (tall canopy: Y = 15.0 to 26.0 mm)
    # Chamber 2: 100 Arabidopsis thaliana seeds (miniature prostrate rosette lawn: Y = 15.0 to 19.5 mm)
    # Leaving open plenum above Chamber 2 (Y = 19.5 to 33.67 mm) directly under Port 2 Central Exhaust
    j_brapa = (y_coords >= 0.0150) & (y_coords <= 0.0260)
    j_arab = (y_coords >= 0.0150) & (y_coords <= 0.0195)
    
    k_canopy_ch1 = (z_coords >= -0.0300) & (z_coords <= -0.0100)
    k_canopy_ch2 = (z_coords >= -0.0070) & (z_coords <= +0.0070)
    k_canopy_ch3 = (z_coords >= +0.0100) & (z_coords <= +0.0300)
    i_canopy = (x_coords >= -0.0100) & (x_coords <= +0.0200)

    canopy_ch1 = (grid == 1) & i_canopy[:, None, None] & j_brapa[None, :, None] & k_canopy_ch1[None, None, :]
    canopy_ch2 = (grid == 1) & i_canopy[:, None, None] & j_arab[None, :, None] & k_canopy_ch2[None, None, :]
    canopy_ch3 = (grid == 1) & i_canopy[:, None, None] & j_brapa[None, :, None] & k_canopy_ch3[None, None, :]

    grid[canopy_ch1 | canopy_ch2 | canopy_ch3] = 6 # PLANT_CANOPY

    print("\n--- Domain Summary ---")
    print(f"Grid shape (Nx, Ny, Nz): {grid.shape}")
    print(f"Voxel pitch: {pitch*1000:.2f} mm")
    print(f"Solid wall voxels (0): {np.sum(grid == 0)}")
    print(f"Air headspace voxels (1): {np.sum(grid == 1)}")
    print(f"Porous substrate voxels (2): {np.sum(grid == 2)}")
    print(f"Inlet port voxels (3): {np.sum(grid == 3)}")
    print(f"Outlet port voxels (4): {np.sum(grid == 4)}")
    print(f"Plant canopy voxels (6): {np.sum(grid == 6)}")
    
    vol_air = (np.sum(grid == 1) + np.sum(grid == 6)) * (pitch * 100)**3
    vol_porous = np.sum(grid == 2) * (pitch * 100)**3
    print(f"Air volume: {vol_air:.2f} cm3")
    print(f"Porous bed volume: {vol_porous:.2f} cm3")
    print(f"Total internal volume: {vol_air + vol_porous:.2f} cm3")

    # Save domain
    os.makedirs('geometry', exist_ok=True)
    np.savez_compressed('geometry/gbc2_grid.npz',
                       grid=grid,
                       pitch=pitch,
                       transform=vox.transform,
                       bounds=vox.bounds)
    print("Saved geometry/gbc2_grid.npz successfully!")
    return grid, pitch, vox.transform

if __name__ == '__main__':
    build_grid(pitch=0.001) # 1 mm pitch
