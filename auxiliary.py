import numpy as np, numpy.linalg as npl
import linecache
from scipy.spatial import Voronoi
from itertools import product
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def get_Wigner_Seitz_BZ(lattice_vectors):
    latt = []
    prefactors = [0., -1., 1.]
    for p in prefactors:
        for u in lattice_vectors:
            latt.append(p * u)
    lattice = []
    for vs in product(latt, latt, latt):
        a = vs[0] + vs[1] + vs[2]
        if not any((a == x).all() for x in lattice):
            lattice.append(a)
    voronoi = Voronoi(lattice)
    bz_facets = []
    bz_ridges = []
    bz_vertices = []
    for pid, rid in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        if(pid[0] == 0 or pid[1] == 0):
            bz_ridges.append(voronoi.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(voronoi.vertices[rid])
            bz_vertices += rid
    bz_vertices = list(set(bz_vertices))
    return voronoi.vertices[bz_vertices], bz_ridges, bz_facets


def _interpolate(start, stop, num_points, endpoint = False):
    if endpoint == False:
        divisor = num_points
    else:
        divisor = num_points - 1
            
    start = np.array(start)
    stop = np.array(stop)
        
    steps = (1.0/divisor) * (stop - start)
    k_points_inter = (steps[:,None] * np.arange(num_points) + start[:,None]).T
        
    return k_points_inter


def kPath(nodes, num_kpoints, lattice = None):
 
    nodes = np.array(nodes)
    num_nodes = nodes.shape[0]
    
    # k_space metric tensor
    if lattice is None:
        k_metric = np.identity(nodes.shape[1])
    else:
        k_metric = np.linalg.inv(np.dot(lattice, lattice.T))

    k_nodes = np.zeros(num_nodes, dtype = float)
    for n in range(num_nodes - 1):
        k_segment = nodes[n+1] - nodes[n]
        k_segment_length = np.sqrt(np.dot(k_segment, np.dot(k_metric,k_segment)))
        k_nodes[n+1] = k_nodes[n] + k_segment_length
    
    # Allocate K points to each segment
    node_index = []
    for n in range(num_nodes - 1):
        ratio_distance = k_nodes[n] / k_nodes[-1]
        node_index.append(int(round(ratio_distance * (num_kpoints - 1))))
    node_index.append(num_kpoints - 1)
        
    # Interpolate K Points to each segment
    k_dist = np.array([])
    k_vec = np.array([])
    for n in range(num_nodes - 1):
        k_dist_segment = np.linspace(k_nodes[n], k_nodes[n+1], node_index[n+1] - node_index[n], endpoint = False)
        k_dist = np.append(k_dist, k_dist_segment)
            
        k_vec_segment = _interpolate(nodes[n], nodes[n+1], node_index[n+1] - node_index[n], endpoint = False)
        k_vec = np.append(k_vec, k_vec_segment)
            
    k_dist = np.append(k_dist, k_nodes[-1])
    k_vec = np.append(k_vec, nodes[-1]).reshape((num_kpoints, nodes.shape[1]))

    return k_vec, k_dist, k_nodes