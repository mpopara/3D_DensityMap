"""Compute the density of a (weighted) MD
traj on a grid, splits the ensemble into
two densities for estiamation of resolution via FSC
"""
import pathlib
import numpy as np
import mdtraj as md
import numba as nb
import mrcfile


@nb.jit(nopython=True)
def compute_density(xyz, weights, dg=0.2, split=True, n_x=100, n_y=100, n_z=100):
    x_min = np.min(xyz[:, :, 0])
    x_max = np.max(xyz[:, :, 0])
    dx = x_max - x_min

    y_min = np.min(xyz[:, :, 1])
    y_max = np.max(xyz[:, :, 1])
    dy = y_max - y_min

    z_min = np.min(xyz[:, :, 2])
    z_max = np.max(xyz[:, :, 2])
    dz = z_max - z_min
    
    if n_x < 0:
        ev = lambda x: int((int(x // 2) + int(x % 2)) * 2)
        n_x = ev(dx // dg + 1)
        n_y = ev(dy // dg + 1)
        n_z = ev(dz // dg + 1)

    g = np.zeros((n_x, n_y, n_z), dtype=np.float32)
    g1 = np.zeros((n_x, n_y, n_z), dtype=np.float32)
    g2 = np.zeros((n_x, n_y, n_z), dtype=np.float32)
    n_frames, n_atoms, n_dim = xyz.shape
    for i in range(n_frames):
        weight = weights[i]
        xs = ((xyz[i, :, 0] - x_min) // dg).astype(np.int32)
        ys = ((xyz[i, :, 1] - y_min) // dg).astype(np.int32)
        zs = ((xyz[i, :, 2] - z_min) // dg).astype(np.int32)
        for j in range(n_atoms):
            g[xs[j], ys[j], zs[j]] += weight
    # Split ensemble in two sub-ensembles for precision estimate via FSC
    # Here sub-ensembles are containing odd and even frames
    for i in range(0, n_frames, 2):
        weight = weights[i]
        xs = ((xyz[i, :, 0] - x_min) // dg).astype(np.int32)
        ys = ((xyz[i, :, 1] - y_min) // dg).astype(np.int32)
        zs = ((xyz[i, :, 2] - z_min) // dg).astype(np.int32)
        for j in range(n_atoms):
            g1[xs[j], ys[j], zs[j]] += weight
    for i in range(1, n_frames, 2):
        weight = weights[i]
        xs = ((xyz[i, :, 0] - x_min) // dg).astype(np.int32)
        ys = ((xyz[i, :, 1] - y_min) // dg).astype(np.int32)
        zs = ((xyz[i, :, 2] - z_min) // dg).astype(np.int32)
        for j in range(n_atoms):
            g2[xs[j], ys[j], zs[j]] += weight
    return g, g1, g2




def write_density(traj, weights, prefix, parent_path, voxel_size=0.2, superpose=False):
    if superpose:
        traj.superpose(traj[0])
        traj.center_coordinates()
    g, h1, h2 = compute_density(traj.xyz, weights, dg=voxel_size)
    with mrcfile.new(parent_path / str(prefix + 'g.mrc'), overwrite=True) as mrc:
         mrc.set_data(g)
         mrc.update_header_from_data()
         mrc.voxel_size = voxel_size * 10.0
    with mrcfile.new(parent_path / str(prefix + 'h1.mrc'), overwrite=True) as mrc:
         mrc.set_data(h1)
         mrc.update_header_from_data()
         mrc.voxel_size = voxel_size * 10.0
    with mrcfile.new(parent_path / str(prefix + 'h2.mrc'), overwrite=True) as mrc:
         mrc.set_data(h2)
         mrc.update_header_from_data()
         mrc.voxel_size = voxel_size * 10.0




################# input data paths and input parameters ######################



topPath = pathlib.Path('C:/User/folder/topology.pdb')
 
trajPath = pathlib.Path('C:/User/folder/trajectory.dcd')
 
weightPath = pathlib.Path('C:/User/folder/conformer_weights.dat')

savePath = pathlib.Path('C:/User/folder/')
            
quantiles_ranges = [(0, 0.25), (0.25, 0.50), (0.5, 0.75), (0.75, 1.00)]



##############################################################################


    
t = md.load(str(trajPath), top=str(topPath))
t.superpose(t[0]) # aligning on first frame in trajectory
# How to align is a choice of user
# Tip: When aligning to external reference structure, in later visualization 
# map and ref.structure will not be aligned, therefore you should translate the reference as follows:
# ref_translated.xyz = ref.xyz[:,:,:]-[x_min,y_min,z_min], where  x_min,y_min,z_min
# are returned from compute_density function  


t.center_coordinates()
weights = np.loadtxt(weightPath).T[1]                 

write_density(t, weights, '3D_DensityMap_', savePath)


# compute 3D density for structures that are in different quantile ranges of cumulative weights

weights_ids = weights.argsort()[::-1]
weights_sorted = weights[weights_ids]
weights_sorted_cum = np.cumsum(weights_sorted)

n_struct = len(weights)
quan_structs = list()
for q in quantiles_ranges:
    q_filter = (weights_sorted_cum >  q[0]) & (weights_sorted_cum <= q[1])
    q_idx = np.where(q_filter)[0]
    quan_structs.append(q_idx)
        
traj_sorted = t[weights_ids]
cum_structs = list()
for i, q in enumerate(quantiles_ranges):
    qts = quan_structs[i]
    fn_qua_prefix = 'posterior_Q%.2f-%.2f_' % (q[0], q[1])
    fn_cum_prefix = 'posterior_C0-%.2f_' % (q[1])
    print("Q: ", q[0], ",", q[1], ": ", len(qts))
   # np.savetxt(savePath / str(fn_qua_prefix + "frames.txt"), weights_ids[qts], fmt='%d')
    write_density(
        traj_sorted[qts],
        weights_sorted[qts], 
        fn_qua_prefix, 
        savePath,
        superpose=True
    )
    cum_structs += qts.tolist()
    write_density(
        traj_sorted[cum_structs],
        weights_sorted[cum_structs], 
        fn_cum_prefix, 
        savePath,
        superpose=True
    )
     
        
        