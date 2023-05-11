import brille
import spglib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
from voronoiBZ import soft_bounded_voronoi, unique_rows, buffer_hull, in_hull
from numpy.linalg import eigvalsh
from fractions import Fraction
from auxiliary import kPath
import sys,ast


def direct2cartesian(point, lattice):
    return np.dot(point, lattice).astype(np.float64)

def cartesian2direct(point, lattice):
    return np.dot(point, np.linalg.inv(lattice)).astype(np.float64)

def decimal2fraction(decimal, max_denominator=100):
    v_decimal2fraction = np.vectorize(lambda x: Fraction(x).limit_denominator(max_denominator))
    return v_decimal2fraction(decimal)

def unique_indices(array_1d):
    sort_indices = np.argsort(array_1d)
    array_1d = np.asarray(array_1d)[sort_indices]
    values, first_indices, *others = np.unique(array_1d, return_index=True)
    indices = np.split(sort_indices, first_indices[1:])
    return values, indices

def generate_dyn(q_direct, reciprocal_lattice, num_atoms, num_cells, C, R, M):
    q_cart = direct2cartesian(q_direct, reciprocal_lattice)

    D_q = np.zeros((num_atoms, 3, num_atoms, 3), dtype=np.complex64)
    for i_atom in range(num_atoms):
        for i_cart in range(3):
            for j_atom in range(num_atoms):
                for j_cart in range(3):
                    for i_cell in range(num_cells):
                        exp_q_dot_r = 0.0
                        for i_im in range(len(R[i_cell, i_atom, j_atom])):
                            exp_q_dot_r += np.exp(-2j * np.pi * q_cart.dot(R[i_cell, i_atom, j_atom][i_im]))
                        exp_q_dot_r = exp_q_dot_r / len(R[i_cell, i_atom, j_atom])
                        D_q[i_atom, i_cart, j_atom, j_cart] += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_q_dot_r / np.sqrt(M[i_atom] * M[j_atom])
    return D_q

def generate_dyn_vornoi(q_direct, reciprocal_lattice, num_atoms, nums_cells, C_list, R_list, M, weights):
    D_q = np.zeros((num_atoms, 3, num_atoms, 3), dtype=np.complex64)
    for i in range(len(nums_cells)):
        D_q += generate_dyn(q_direct, reciprocal_lattice, num_atoms, nums_cells[i], C_list[i], R_list[i], M) * weights[i]
    return D_q

def get_savt_weight(points, voronoi_cells, ibz_point_indicators):
    savt_weight = np.zeros(max(ibz_point_indicators) + 1)
    for i, point in enumerate(points):
        if tuple(point) in voronoi_cells:
            savt_weight[ibz_point_indicators[i]] += voronoi_cells[tuple(point)].volume
        else:
            print("Function get_savt_weight cannot find point " + str(point) + " in the given voronoi_cells")

    savt_weight = np.array(savt_weight) / sum(savt_weight)
    return savt_weight

def assign_points2grids(points_direct, grids_sizes):
    points_direct = decimal2fraction(points_direct)
    grid_indicators = np.full(len(points_direct), np.nan)

    flag = 0
    for grids_size in grids_sizes:
        for i, point in enumerate(points_direct):
            denominators = [coordinate.denominator for coordinate in point]
            if (np.gcd(denominators, grids_size) == denominators).all():
                grid_indicators[i] = flag
        flag = flag + 1
    return grid_indicators.astype(np.int64)

def mark_symmetry_images(points_con, brille_bz):
    if points_con.shape[1] == 2:
        padding = decimal2fraction(np.zeros((points_con.shape[0], 1)))
        points_con_3d = np.hstack((points_con, padding))
    else:
        points_con_3d = points_con

    ref = decimal2fraction(unique_rows(brille_bz.ir_moveinto(points_con_3d)[0]))
    ibz_point_indicators = np.full(len(points_con_3d), np.nan)

    points_ibz = decimal2fraction(brille_bz.ir_moveinto(points_con_3d)[0])
    
    for i, point in enumerate(points_ibz):
        ibz_point_indicators[i] = np.where(np.all(point==ref, axis=1))[0][0]

    return ibz_point_indicators.astype(np.int64)

def extend_bz_points(points_bz_direct):
    dim = points_bz_direct.shape[1]

    if dim == 2:
        points_direct = np.array([]).reshape(0, 2)
        for i in [0, 1, -1]:
            for j in [0, 1, -1]:
                points_direct = np.vstack((points_direct, points_bz_direct + np.array([i, j])))

    else:
        points_direct = np.array([]).reshape(0, 3)
        for i in [0, 1, -1]:
            for j in [0, 1, -1]:
                for k in [0, 1, -1]:
                    points_direct = np.vstack((points_direct, points_bz_direct + np.array([i, j, k])))
    return points_direct

dim = 3
lte_list = ast.literal_eval(sys.argv[1])
lte_list = [str(i) for i in lte_list]

# print(lte_list,'lte')
path = './'
num_atoms = np.loadtxt(lte_list[0]+"/equilibrium.dat", dtype=np.int64, comments=['#', '$', '@'], max_rows=1)
M = np.loadtxt(lte_list[0]+"/equilibrium.dat", dtype=np.float64, comments=['#', '$', '@'], skiprows=1, usecols=1)

#Regardless of the dimensionality of the system, one needs a 3D lattice so that spglib can work
lattice_given_3d = np.loadtxt(lte_list[0]+"/lattice.dat", dtype=np.float64, comments=['#', '$', '@'])#, usecols=range(dim), max_rows=dim)
positions_given_cart = np.loadtxt(lte_list[0]+"/equilibrium.dat", dtype=np.float64, comments=['#', '$', '@'], skiprows=1, usecols=(2,3,4))
atom_types_given_str = np.loadtxt(lte_list[0]+"/equilibrium.dat", dtype=str, comments=['#', '$', '@'], skiprows=1, usecols=0)

positions_given_direct = cartesian2direct(positions_given_cart, lattice_given_3d)
type2int = dict([(y, x) for x, y in enumerate(sorted(set(atom_types_given_str)))])
atom_types_given_int = [type2int[i] for i in atom_types_given_str]

cell_given = (lattice_given_3d, positions_given_direct, atom_types_given_int)
symmetry_data = spglib.get_symmetry_dataset(cell_given)

# spglib gives us the standard primitive and conventional lattice
lattice_pri, _, _ = spglib.standardize_cell(cell_given, to_primitive=True, no_idealize=False, symprec=1e-5)
lattice_con, _, _ = spglib.standardize_cell(cell_given, to_primitive=False, no_idealize=False, symprec=1e-5)

# obtain the transformation matrix from the primitive lattice to the conventional lattice
pri2con = decimal2fraction(np.dot(np.linalg.inv(lattice_pri), lattice_con))
con2pri = decimal2fraction(np.dot(np.linalg.inv(lattice_con), lattice_pri))

symmetry_indicator = symmetry_data["international"]
print("Spglib finds the space group: " + str(symmetry_indicator))
point_group = symmetry_data["rotations"]

# create Brille lattice object with correct symmetry. Remember for automatical IBZ searching Brille requires the conventional lattice instead of the primitive
lattice_brille = brille.Direct(lattice_con, symmetry_indicator)

# find the BZ and the IBZ (automatically)
reciprocal_lattice_brille = lattice_brille.star
bz = brille.BrillouinZone(reciprocal_lattice_brille, time_reversal_symmetry=True)


grids_sizes = []
q_points_3d = []
C_list = []
R_list = []

for indicator in lte_list:
    grids_size = np.loadtxt(path+indicator+"/grid.dat", dtype=np.int64, comments=['#', '$', '@'], usecols=range(dim), max_rows=1)
    grids_sizes.append(grids_size)

    #Regardless of the dimensionality of the system, one needs a 3D q points so that brille can work
    q_points_per_grid = np.loadtxt(path+indicator+"/ibz.dat", dtype=np.float64, comments=['#', '$', '@'], usecols=range(3))
    q_points_3d.append(q_points_per_grid)

    C = np.zeros((np.prod(grids_size), num_atoms, 3, num_atoms, 3), dtype=np.complex64)
    with open(path+indicator+"/force.dat",'r') as f:
        for line in f:
            C[int(line.split()[0])-1, int(line.split()[1])-1, int(line.split()[2])-1, int(line.split()[3])-1, int(line.split()[4])-1] = -float(line.split()[5])
    C_list.append(C)
    f.close()

    R = np.frompyfunc(list, 0, 1)(np.empty((np.prod(grids_size), num_atoms, num_atoms), dtype=object))
    with open(path+indicator+"/delta_prim.dat",'r') as f:
        if dim == 2:
            for line in f:
                R[int(line.split()[0])-1, int(line.split()[1])-1, int(line.split()[2])-1].append(np.array([float(line.split()[4]), float(line.split()[5])]))
        else:
            for line in f:
                R[int(line.split()[0])-1, int(line.split()[1])-1, int(line.split()[2])-1].append(np.array([float(line.split()[4]), float(line.split()[5]), float(line.split()[6])]))
    R_list.append(R)
    f.close()

q_points_3d = unique_rows(np.vstack((q_points_3d)))
nums_cells = list(map(np.prod, grids_sizes))

print("Finish reading all lte folders")



# exit()


# Convert the coordinate of the input q points from the primitive lattice (spglib used) to the conventional lattice (Brille used)
q_points_con = np.dot(q_points_3d, pri2con)
q_points_ibz_con = decimal2fraction(bz.ir_moveinto(q_points_con)[0])
q_points_ibz_pri = np.dot(q_points_ibz_con, con2pri)

q_points_bz_pri = np.array([]).reshape(0, 3)

for rotation in point_group:
    q_points_bz_pri = np.vstack((q_points_bz_pri, np.dot(q_points_ibz_pri, rotation)))

q_points_bz_pri = unique_rows(q_points_bz_pri)

print("q_points_ibz_pri and q_points_bz_pri are ready")




# For Voronoi tessellation, one needs to use the correct dimensionality (2D or 3D)
lattice = lattice_pri[0:dim, 0:dim]
reciprocal_lattice = np.linalg.inv(lattice).T

boundary_con = unique_rows(bz.ir_vertices[:, 0:dim])
boundary_pri = np.dot(boundary_con, con2pri[0:dim, 0:dim])
boundary = direct2cartesian(boundary_pri, reciprocal_lattice)

q_points_ibz_direct = unique_rows(q_points_ibz_pri[:, 0:dim])
q_points_ibz = direct2cartesian(q_points_ibz_direct, reciprocal_lattice)

q_points_bz_direct = unique_rows(q_points_bz_pri[:, 0:dim] % 1)
q_points_bz = direct2cartesian(q_points_bz_direct, reciprocal_lattice)

q_points_direct = extend_bz_points(q_points_bz_direct)
q_points = direct2cartesian(q_points_direct, reciprocal_lattice)
print("Find " + str(len(q_points)) + " q points in the extended BZ")

buffer = buffer_hull(boundary, r=None, num_pts=10) #Boundaries are NOT mirror planes
q_points = q_points[in_hull(q_points, buffer)]
q_points_direct = decimal2fraction(cartesian2direct(q_points, reciprocal_lattice))
print("Consider " + str(len(q_points)) + " q points in the buffered BZ")

grid_indicators = assign_points2grids(q_points_direct, grids_sizes)
ibz_point_indicators = mark_symmetry_images(np.dot(q_points_direct, pri2con[0:dim, 0:dim]), bz)

voronoi_cells = soft_bounded_voronoi(q_points, boundary, r=None)
q_points_savt = np.array(list(voronoi_cells.keys()))
q_points_savt_direct = decimal2fraction(cartesian2direct(q_points_savt, reciprocal_lattice))
savt_point_indicators = mark_symmetry_images(np.dot(q_points_savt_direct, pri2con[0:dim, 0:dim]), bz)
savt_grid_indicators = assign_points2grids(q_points_savt_direct, grids_sizes)
savt_weight = get_savt_weight(q_points_savt, voronoi_cells, savt_point_indicators)
print("There are " + str(len(q_points_savt)) + " q points with none zero SAVT weights")
print(savt_weight / savt_weight[0])

# fig = plt.figure()
# ax = Axes3D(fig)
# ax._axis3don = False

# cmap = cm.brg
# norm = Normalize(vmin=min(savt_grid_indicators), vmax=max(savt_grid_indicators))

# for i in range(len(lte_list)):
#     ax.scatter(q_points_savt[savt_grid_indicators==i, 0], q_points_savt[savt_grid_indicators==i, 1], q_points_savt[savt_grid_indicators==i, 2], marker='x', color=cmap(norm(i)))

# i = 0
# # print(voronoi_cells.items())
# # for point, voronoi_cell in voronoi_cells.items():
# #     for s in voronoi_cell.simplices:
# #         print('SSS',s,np.shape(voronoi_cell.points[s]))
# #         cell = Poly3DCollection(voronoi_cell.points[s])
# #         cell.set_facecolor(cmap(norm(savt_grid_indicators[i])))
# #         cell.set_alpha(0.15)
# #         ax.add_collection3d(cell)
# #     i = i + 1

# # plt.show()


# exit()

grid_target = sys.argv[2]

q_target= np.loadtxt(grid_target + "/ibz.dat", dtype=np.float64, comments=['#', '$', '@'], usecols=range(3))
q_target_con = np.dot(q_target, pri2con)
q_target_ibz_con = decimal2fraction(bz.ir_moveinto(q_target_con)[0])
q_target_ibz_pri = np.dot(q_target_ibz_con, con2pri)

q_target_bz_pri = np.array([]).reshape(0, 3)

for rotation in point_group:
    q_target_bz_pri = np.vstack((q_target_bz_pri, np.dot(q_target_ibz_pri, rotation)))

q_target_bz_pri = unique_rows(q_target_bz_pri)
q_target_bz_con = np.dot(q_target_bz_pri, pri2con)

print("finish reading the target lte folder")



#Calculate dynamical matrix
values, indices = unique_indices(grid_indicators)
print('values',values, indices)
for i, q_interp_ibz_pri in enumerate(q_target_ibz_pri):
    q_interp_ibz_direct = q_interp_ibz_pri[0:dim]
    q_interp_ibz = direct2cartesian(q_interp_ibz_direct, reciprocal_lattice)
    points_weight = np.zeros(len(q_points))
    
    overlap = np.where(np.all(np.isclose(q_interp_ibz, q_points), axis=1))[0]

    if overlap.size == 0:
        voronoi_cells_ref = soft_bounded_voronoi(q_points, boundary, r=None)
    
        q_interp_bz_pri = np.array([]).reshape(0, 3)
        for rotation in point_group:
            q_interp_bz_pri = np.vstack((q_interp_bz_pri, np.dot(q_interp_ibz_pri, rotation)))
        q_interp_bz_direct = unique_rows(q_interp_bz_pri[:, 0:dim] % 1)

        q_interp_direct = extend_bz_points(q_interp_bz_direct)
        q_interp_direct = q_interp_direct[in_hull(direct2cartesian(q_interp_direct, reciprocal_lattice), buffer)]
        
        q_interp_with_images = direct2cartesian(q_interp_direct, reciprocal_lattice)

        points_with_interp = np.vstack((q_points, q_interp_with_images))
        voronoi_cells_interp = soft_bounded_voronoi(points_with_interp, boundary, r=None) #Boundaries are NOT mirror planes

        #double check here
        for point, voronoi_cell in voronoi_cells_ref.items():
            index = np.where(np.all(np.isclose(point, q_points), axis=1))[0][0]
            print('index',index)
            print('volume',voronoi_cell.volume)
            print('point',point,len(voronoi_cells_interp))
            print(voronoi_cells_interp[point].volume)
            points_weight[index] = (voronoi_cell.volume - voronoi_cells_interp[point].volume)

        points_weight = points_weight / sum(points_weight)

    else:
        points_weight[overlap[0]] = 1.0

    grids_weight = np.zeros(len(lte_list))
    for value, index in zip(values, indices):
        grids_weight[value] = points_weight[index].sum()

    D_q = generate_dyn_vornoi(q_interp_ibz_direct, reciprocal_lattice, num_atoms, nums_cells, C_list, R_list, M, grids_weight)
    print(q_interp_ibz_direct)
    print(grids_weight)
    with open(grid_target +  "/dyn_mat." + str(i+1) + ".dat", 'a') as f:
        for i_atom in range(num_atoms):
            for i_cart in range(3):
                for j_atom in range(num_atoms):
                    for j_cart in range(3):
                        DD=D_q[i_atom, i_cart, j_atom, j_cart]
                        f.write("%d %d %d %d %.17g %.17g\n" %(i_atom+1, i_cart+1, j_atom+1 ,j_cart+1, -np.real(DD), -np.imag(DD)))
    f.close()

    with open(grid_target +  "/atoms_in_primitive_cell." + str(i+1) + ".dat", 'a') as f:
        for i in range(num_atoms):
            f.write("%s %.17g %.17g %.17g %.17g\n" %(atom_types_given_str[i], M[i], positions_given_direct[i, 0], positions_given_direct[i, 1], positions_given_direct[i, 2]))
    f.close()

print("Calculating weights for grid " + grid_target)
q_target_bz_pri = np.array([]).reshape(0, 3)
for rotation in point_group:
    q_target_bz_pri = np.vstack((q_target_bz_pri, np.dot(q_target_ibz_pri, rotation)))

q_target_bz_direct = unique_rows(q_target_bz_pri[:, 0:dim] % 1)
q_target_direct = extend_bz_points(q_target_bz_direct)

q_target = direct2cartesian(q_target_direct, reciprocal_lattice)

voronoi_cells = soft_bounded_voronoi(q_target, boundary, r=None) #Boundaries are NOT mirror planes
q_points = np.array(list(voronoi_cells.keys()))

q_points_direct = decimal2fraction(cartesian2direct(q_points, reciprocal_lattice))
grid_indicators = assign_points2grids(q_points_direct, grids_sizes)
ibz_point_indicators = mark_symmetry_images(np.dot(q_points_direct, pri2con[0:dim, 0:dim]), bz)

savt_weight = get_savt_weight(q_points, voronoi_cells, ibz_point_indicators)
print(savt_weight/savt_weight[0])

with open(grid_target + "/lte" + "/ibz.dat", 'w') as f:
    for i in range(len(q_target_ibz_pri[:, 0:dim])):
        multiplicity = round(savt_weight[i] / savt_weight[0])
        f.write("%.17g %.17g %.17g %d\n" %(q_target_ibz_pri[i, 0], q_target_ibz_pri[i, 1], q_target_ibz_pri[i, 2], multiplicity))
f.close()

with open(grid_target + "/lte" + "/kpoint_to_supercell.dat", 'a') as f:
    for i in range(len(q_target_ibz_pri[:, 0:dim])):
        f.write("%.17g %.17g %.17g %d\n" %(q_target_ibz_pri[i, 0], q_target_ibz_pri[i, 1], q_target_ibz_pri[i, 2], i+1))
f.close()

print('All Job Done')