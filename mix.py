from farey2qe import *
import sys, ast

def get_D(script_dir):
    #get the dynamical matrix at gamma point for smaller cells
    num_atoms = np.loadtxt(lte_list[0] + "/equilibrium.dat", dtype=np.int64, comments=['#', '$', '@'], max_rows=1)
    data = np.loadtxt(script_dir + '/force.dat')
    indices = data[:, :5].astype(int) - 1  # subtract 1 to convert to 0-based indexing
    values = data[:, 5]
    # determine shape of 5D array
    shape = np.max(indices, axis=0) + 1
    # create 5D array and fill with values
    C = np.empty(shape)
    C[indices[:, 0], indices[:, 1], indices[:, 2], indices[:, 3], indices[:, 4]] = -values

    data = np.loadtxt(script_dir + '/delta_prim.dat')
    indices = data[:, :4].astype(int) - 1  # subtract 1 to convert to 0-based indexing
    values = data[:, 4:]
    shape = np.max(indices, axis=0) + 1
    R = np.empty(shape, dtype=list)
    R[indices[:, 0], indices[:, 1], indices[:, 2], indices[:, 3]] = values.tolist()
    q=np.array([0,0,0])
    D = generate_dyn_qe(q , C, R, 1)
    return  np.reshape(D, (-1, 3*num_atoms))

def mix(D3,D4,alpha):
    
    header='''Dynamical matrix file
File generated with the CellConstructor by Lorenzo Monacelli
2 2 0     8.7523599000000054     0.0000000000000000     0.0000000000000000     0.0000000000000000     0.0000000000000000     0.0000000000000000
Basis vectors
   -0.0000000000000000     0.6898687370000000     0.6898687370000000
    0.6898687370000000    -0.0000000000000000     0.6898687370000000
    0.6898687370000000     0.6898687370000000    -0.0000000000000000
        1  'Sn '  108197.5460994285967899
        2  'Te '  116300.2854206645570230
    1     1    -0.0000000000000000     0.0000000000000000    -0.0000000000000000
    2     2     0.6898687368000000     0.6898687368000000     0.6898687368000000
'''

    with open(f'./{lte_list[-1]}/{dyn_name[-1]}1','w') as f:
        f.write(header)
        D = alpha*D3 + (1-alpha)*D4
        D_real = np.real(D); D_imag = np.imag(D)
        f.write('\n     Dynamical Matrix in cartesian axes\n')
        f.write('\n')
        f.write(f'    q = (     0.000000000000     0.000000000000     0.000000000000 )\n')
        f.write('\n')
        
        for j in range(num_atoms):
            for k in range(num_atoms):
                f.write(f'    {j+1}    {k+1}\n')
                for l in range(3):
                    print(j,k,l)
                    f.write(f'  {D_real[3*j+l,3*k]:>14.9f}{D_imag[3*j+l,3*k]:>14.9f}{D_real[3*j+l,3*k+1]:>14.9f}{D_imag[3*j+l,3*k+1]:>14.9f}{D_real[3*j+l,3*k+2]:>14.9f}{D_imag[3*j+l,3*k+2]:>14.9f}\n')
        charge=get_charge(f'./{lte_list[0]}/{dyn_name[0]}',1)
        f.write('\n')
        if charge:
            for ll in charge:
                f.write(f'{ll}')

if __name__=='__main__':
    lte_list = ast.literal_eval(sys.argv[1])
    lte_list = [str(i) for i in lte_list]
    dyn_name = ast.literal_eval(sys.argv[2])
    dyn_name = [str(i) for i in dyn_name]
    alat=1
    D3=get_D(lte_list[0])
    D4=get_D(lte_list[1])
    alpha=float(sys.argv[3])
    num_atoms = np.loadtxt(lte_list[0] + "/equilibrium.dat", dtype=np.int64, comments=['#', '$', '@'], max_rows=1)
    mix(D3,D4,alpha)