
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from auxiliary import kPath



def generate_dyn(k, reciprocal_lattice, num_cell, num_atoms, C, R, M):

    k_abs = np.dot(k, np.transpose(reciprocal_lattice))
    # print(k_abs)
    D_q = np.zeros((num_atoms, 3, num_atoms, 3), dtype=complex)
    for i_atom in range(num_atoms):
        for i_cart in range(3):
            for j_atom in range(num_atoms):
                for j_cart in range(3):
                    for i_cell in range(num_cell):
                        exp_k_dot_r = 0.0
                        for i_im in range(len(R[i_cell, i_atom, j_atom])):
                            exp_k_dot_r += np.exp(-2j*np.pi*k_abs.dot(R[i_cell,i_atom,j_atom][i_im]))
                        exp_k_dot_r = exp_k_dot_r / len(R[i_cell, i_atom, j_atom])
                        D_q[i_atom, i_cart, j_atom, j_cart] += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_k_dot_r /np.sqrt(M[i_atom]*M[j_atom])
    return D_q

def debug_x(k, reciprocal_lattice, num_cell,i_atom,j_atom, C, R, M):

    k_abs = np.dot(k, np.transpose(reciprocal_lattice))
    D_q = 0

    for i_cart in [0]:
                for j_cart in [0]:
                    for i_cell in range(num_cell):
                        exp_k_dot_r = 0.0
                        for i_im in range(len(R[i_cell, i_atom, j_atom])):
                            exp_k_dot_r += np.exp(-2j*np.pi*k_abs.dot(R[i_cell,i_atom,j_atom][i_im]))
                        exp_k_dot_r = exp_k_dot_r / len(R[i_cell, i_atom, j_atom])
                        # exp_k_dot_r = 1
                        D_q += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_k_dot_r /np.sqrt(M[i_atom]*M[j_atom])
                        print(i_cell,C[i_cell, i_atom, i_cart, j_atom, j_cart],exp_k_dot_r,C[i_cell, i_atom, i_cart, j_atom, j_cart]* exp_k_dot_r)
    return D_q

def debug_y(k, reciprocal_lattice, num_cell,i_atom,j_atom, C, R, M):

    k_abs = np.dot(k, np.transpose(reciprocal_lattice))
    D_q = 0

    for i_cart in [1]:
                for j_cart in [1]:
                    for i_cell in range(num_cell):
                        exp_k_dot_r = 0.0
                        for i_im in range(len(R[i_cell, i_atom, j_atom])):
                            exp_k_dot_r += np.exp(-2j*np.pi*k_abs.dot(R[i_cell,i_atom,j_atom][i_im]))
                        exp_k_dot_r = exp_k_dot_r / len(R[i_cell, i_atom, j_atom])
                        # exp_k_dot_r = 1
                        D_q += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_k_dot_r /np.sqrt(M[i_atom]*M[j_atom])
                        print(i_cell,C[i_cell, i_atom, i_cart, j_atom, j_cart],exp_k_dot_r,C[i_cell, i_atom, i_cart, j_atom, j_cart]* exp_k_dot_r)
    return D_q

def get_q(file):
    q=[]

    with open(file,'r') as f:
        l=f.readlines()
        for i in range(len(l)):
            if 'Basis vectors' in l[i]:
                L = np.loadtxt(l[i+1:i+4],dtype=float)
                break
        for i in range(len(l)):
            if l[i].strip().startswith('q ='):
                 q.append(np.dot(np.array(l[i].split()[3:6],dtype=float),np.transpose(L)))
    return np.array(q,dtype=float)[:-1]


def write_dynamical():
    ibz=np.loadtxt('ibz.dat')
    nq=len(ibz)
    nq=1
    for i in range(1,nq+1):
        with open(f'dyn_666_{i}','w') as f:
            f.write('''Dynamical matrix file

  2    2   0   8.7523599   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
Basis vectors
      0.000000000    0.709233838    0.709233838
      0.709233838    0.000000000    0.709233838
      0.709233838    0.709233838   -0.000000000
           1  'Pb  '    188851.24717211360
           2  'Te  '    116300.28542066456
    1    1      0.0000000000     -0.0000000000     -0.0000000000
    2    2      0.7092338381      0.7092338381      0.7092338381
            
              ''')
            q = get_q(f'./q-point/harmonic_222_dyn{i}')
            for j in range(len(q)):
                f.write('\nDynamical  Matrix in cartesian axes\n')
                f.write('\n')
                f.write(f'   q = ( {q[j,0]} {q[j,1]} {q[j,2]} )\n')
                D_q = generate_dyn(q[j], reciprocal_lattice, num_cell, num_atoms, C, R, M)
                D = np.reshape(D_q, (num_atoms*3, num_atoms*3))
                D = (np.matrix(D) + np.matrix(D).H)/2
                print(np.shape(D))
                print(D)
                D_real = np.real(D); D_imag = np.imag(D)
                for jj in range(2):
                    for kk in range(2):
                        f.write(f'    {jj+1}    {kk+1}\n')
                        for l in range(3):
                            # print('a',D_real[3*jj+l,3*kk])
                            # print('b',D_imag[3*jj+l,3*kk])
                            f.write(f'  {D_real[3*jj+l,3*kk]:>20.15f}{D_imag[3*jj+l,3*kk]:>20.15f}{D_real[3*jj+l,3*kk+1]:>20.15f}{D_imag[3*jj+l,3*kk+1]:>20.15f}{D_real[3*jj+l,3*kk+2]:>20.15f}{D_imag[3*jj+l,3*kk+2]:>20.15f}\n')









if __name__=="__main__":
    script_dir = sys.argv[1]
    dim = 3
    lattice = np.loadtxt(script_dir + "/lattice.dat", dtype=np.float64, comments=['#', '$', '@'])#, usecols=range(dim), max_rows=dim)
    reciprocal_lattice = np.linalg.inv(lattice).T

    num_atoms = np.loadtxt(script_dir + "/equilibrium.dat", dtype=np.int64, comments=['#', '$', '@'], max_rows=1)
    M = np.loadtxt(script_dir + "/equilibrium.dat", dtype=np.float64, comments=['#', '$', '@'], skiprows=1, usecols=1)
    # M = np.array([1,1])
    grids_size = np.loadtxt(script_dir + "/grid.dat", dtype=np.int64, comments=['#', '$', '@'], usecols=range(dim), max_rows=1)
    num_cell = np.prod(grids_size)

    C=np.zeros((num_cell, num_atoms, 3, num_atoms, 3), dtype=complex)

    with open(script_dir + "/force.dat",'r') as f:
        for line in f:
            C[int(line.split()[0])-1, int(line.split()[1])-1, int(line.split()[2])-1, int(line.split()[3])-1, int(line.split()[4])-1] = -float(line.split()[5])

    R=np.frompyfunc(list, 0, 1)(np.empty((num_cell, num_atoms, num_atoms), dtype=object))
    with open(script_dir + "/delta_prim.dat",'r') as f:
        for line in f:
            R[int(line.split()[0])-1, int(line.split()[1])-1, int(line.split()[2])-1].append(np.array([float(line.split()[4]), float(line.split()[5]), float(line.split()[6])]))


    # band structure
    high_sym_points = [[0,0.0,0], [0, 0, 0.5],  [0, .5, .5], [0,0,0], [0.5, 0.5, 0.5]]
    num_kpoints = 200
    k_path, k_distance, k_nodes = kPath(high_sym_points, num_kpoints, lattice)


    # k_path=[[0,0.5,0.5],[0.5,0,0.5]]
    # for k in k_path:
    #     D_q = generate_dyn(k, reciprocal_lattice, num_cell, num_atoms, C, R, M)
    #     D = np.reshape(D_q, (num_atoms*3, num_atoms*3))
    #     D = (np.matrix(D) + np.matrix(D).H)/2
        # print(D)
    # print(get_q('./q-point/harmonic_222_dyn2'))
    # write_dynamical()
    # exit()



    eig_val = []
    for k in k_path:
        D_q = generate_dyn(k, reciprocal_lattice, num_cell, num_atoms, C, R, M)
        D = np.reshape(D_q, (num_atoms*3, num_atoms*3))
        D = (np.matrix(D) + np.matrix(D).H)/2
        w2, v = np.linalg.eigh(D)
        w = np.real(np.sqrt(w2.astype(dtype=complex))) - np.imag(np.sqrt(w2.astype(dtype=complex)))
        print(w);eig_val.append(w)




    eig_val = np.array(eig_val)
    fig, ax = plt.subplots()
    ax.set_xlim(k_nodes[0], k_nodes[-1])
    ax.set_xticks(k_nodes)
    label = ( r'$\Gamma $', r'$x$', r'$M$', r'$\Gamma $', r'$L$')
    ax.set_xticklabels(label)
    for n in range(len(k_nodes)):
        ax.axvline(x = k_nodes[n], linewidth = 0.8, color='k')
    for band_index in range(eig_val.shape[1]):
        ax.plot(k_distance, eig_val[:, band_index], linewidth = 2)

    plt.show()

