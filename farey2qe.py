# from testfourier import *
import numpy as np
import sys

def generate_dyn_qe(k, C, R, alat):
	#Note that k is in cartesian coordinates and needs to be devided by alat
	num_cell, num_atoms, _, _, _ = C.shape
	k_cart = k / alat

	D_q = np.zeros((num_atoms, 3, num_atoms, 3), dtype=complex)
	for i_atom in range(num_atoms):
		for i_cart in range(3):
			for j_atom in range(num_atoms):
				for j_cart in range(3):
					for i_cell in range(num_cell):
						exp_k_dot_r = 0.0
						r = R[i_cell,i_atom,j_atom,:]
						num_image = np.count_nonzero(r != None)
						for i_im in range(num_image):
							# print('kdot',k_cart.dot(r[i_im]))
							exp_k_dot_r += np.exp(-2j * np.pi * k_cart.dot(r[i_im]))
						exp_k_dot_r = exp_k_dot_r / num_image
						D_q[i_atom, i_cart, j_atom, j_cart] += C[i_cell, i_atom, i_cart, j_atom, j_cart] * exp_k_dot_r
	# tricky index switching that is related to the definition of qe matrix
	D_q = np.transpose(D_q, (2, 3, 0, 1)) * 2  # Hatree to Rydberg
	return D_q



def qpoint(input_path,index):
	with open(f'{input_path}{index}','r') as f:
		lines = f.readlines()
	q=[]
	for i in range(len(lines)):
		if ('Dynamical Matrix in cartesian axes' in lines[i]) or ('Dynamical  Matrix in cartesian axes' in lines[i]):
			q.append(lines[i+2].split()[3:6])
			
	return np.array(q,dtype=float)

def get_charge(i_path,index):
	with open(f'{i_path}{index}','r') as f:
		lines = f.readlines()
	for i in range(len(lines)):
		
		if 'Dielectric Tensor' in lines[i]:
			
			return lines[i:i+17]

def get_reduced_lattice():
	with open(f'./{path}/harmonic_{supercell}_dyn1','r') as f3:
		lines=f3.readlines()
		for i in range(len(lines)):
			if 'Basis vectors' in lines[i]:
				L = np.loadtxt(lines[i+1:i+4])
				return L

def output(q,ipath,opath,index):
	global reciprocal_lattice, num_cell, num_atoms, C, R, M
	header='''Dynamical matrix file
!ScVSn
  3   13   0  10.3310004   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
Basis vectors
      1.000000000    0.000000000    0.000000000
     -0.500000000    0.866025404    0.000000000
      0.000000000    0.000000000    1.675419294
           1  'Sc  '    40974.8071860994
           2  'V   '    46430.3369103196
           3  'Sn  '    108197.546099429
    1    1      0.0000000000      0.0000000000      0.0000000000
    2    2      0.5000000000      0.0000000000      0.4156210873
    3    2     -0.2500000000      0.4330127019      0.4156210873
    4    2      0.2500000000      0.4330127019      0.4156210873
    5    2     -0.2500000000      0.4330127019      1.2597982066
    6    2      0.5000000000      0.0000000000      1.2597982066
    7    2      0.2500000000      0.4330127019      1.2597982066
    8    3      0.0000000000      0.0000000000      1.1407148068
    9    3      0.0000000000      0.0000000000      0.5347044871
   10    3      0.5000000000      0.2886751346      0.0000000000
   11    3     -0.0000000001      0.5773502692      0.0000000000
   12    3      0.5000000000      0.2886751346      0.8377096469
   13    3     -0.0000000001      0.5773502692      0.8377096469
'''

	with open(f'{opath}{index}','w') as f:
		f.write(header)
		for i in range(len(q)):
			# D = generate_dyn_qe_1(q[i] , C)
			D = generate_dyn_qe(q[i] , C, R, alat)
			D = np.reshape(D, (-1, 3*num_atoms))
			D = (np.matrix(D) + np.matrix(D).H)/2
			w2, v = np.linalg.eigh(D)
			w = np.real(np.sqrt(w2.astype(dtype=complex))) - np.imag(np.sqrt(w2.astype(dtype=complex)))
			print(w)
			D_real = np.real(D); D_imag = np.imag(D)
			f.write('\n     Dynamical  Matrix in cartesian axes\n')
			f.write('\n')
			f.write(f'     q = ( {q[i][0]} {q[i][1]} {q[i][2]} )\n')
			f.write('\n')
			for j in range(num_atoms):
				for k in range(num_atoms):
					f.write(f'    {j+1}    {k+1}\n')
					for l in range(3):
						f.write(f'  {D_real[3*j+l,3*k]:>14.9f}{D_imag[3*j+l,3*k]:>14.9f}{D_real[3*j+l,3*k+1]:>14.9f}{D_imag[3*j+l,3*k+1]:>14.9f}{D_real[3*j+l,3*k+2]:>14.9f}{D_imag[3*j+l,3*k+2]:>14.9f}\n')
		charge=get_charge(ipath,index)
		f.write('\n')
		if charge:
			for ll in charge:
				f.write(f'{ll}')

def output_new(q,ipath,opath,index):
	global reciprocal_lattice, num_cell, num_atoms, C, R, M
	header='''Dynamical matrix file
!ScVSn
  3   13   0  10.3310004   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000
Basis vectors
      1.000000000    0.000000000    0.000000000
     -0.500000000    0.866025404    0.000000000
      0.000000000    0.000000000    1.675419294
           1  'Sc  '    40974.8071860994
           2  'V   '    46430.3369103196
           3  'Sn  '    108197.546099429
    1    1      0.0000000000      0.0000000000      0.0000000000
    2    2      0.5000000000      0.0000000000      0.4156210873
    3    2     -0.2500000000      0.4330127019      0.4156210873
    4    2      0.2500000000      0.4330127019      0.4156210873
    5    2     -0.2500000000      0.4330127019      1.2597982066
    6    2      0.5000000000      0.0000000000      1.2597982066
    7    2      0.2500000000      0.4330127019      1.2597982066
    8    3      0.0000000000      0.0000000000      1.1407148068
    9    3      0.0000000000      0.0000000000      0.5347044871
   10    3      0.5000000000      0.2886751346      0.0000000000
   11    3     -0.0000000001      0.5773502692      0.0000000000
   12    3      0.5000000000      0.2886751346      0.8377096469
   13    3     -0.0000000001      0.5773502692      0.8377096469
'''

	with open(f'{opath}{index}','w') as f:
		f.write(header)
		for i in range(len(q)):
			# D = generate_dyn_qe_1(q[i] , C)
			D = generate_dyn_qe(q[i] , C, R, alat)
			D_real = np.real(D); D_imag = np.imag(D)
			f.write('\n     Dynamical  Matrix in cartesian axes\n')
			f.write('\n')
			f.write(f'     q = ( {q[i][0]} {q[i][1]} {q[i][2]} )\n')
			f.write('\n')
			for i_atom in range(num_atoms):
				for j_atom in range(num_atoms):
					f.write(f'    {i_atom+1}    {j_atom+1}\n')
					for i_cart in range(3):
						# f.write(f'  {D_real[j,l,k,0]:>14.9f}{D_imag[j,l,k,0]:>14.9f}{D_real[j,l,k,1]:>14.9f}{D_imag[j,l,k,1]:>14.9f}{D_real[j,l,k,2]:>14.9f}{D_imag[j,l,k,2]:>14.9f}\n')
						for j_cart in range(3):
							f.write(f'  {D_real[i_atom,i_cart,j_atom,j_cart]:>14.9f}{D_imag[i_atom,i_cart,j_atom,j_cart]:>14.9f}')
						f.write('\n')
		charge=get_charge(ipath,index)
		f.write('\n')
		if charge:
			for ll in charge:
				f.write(f'{ll}')

if __name__=='__main__':
	
	script_dir = sys.argv[1]; path = script_dir
	supercell = sys.argv[1].split('_')[-1]
	dim = 3		# dimension of the system	
	lattice = np.loadtxt(script_dir + "/lattice.dat", dtype=np.float64, comments=['#', '$', '@'])#, usecols=range(dim), max_rows=dim)

	num_atoms = np.loadtxt(script_dir + "/equilibrium.dat", dtype=np.int64, comments=['#', '$', '@'], max_rows=1)
	M = np.loadtxt(script_dir + "/equilibrium.dat", dtype=np.float64, comments=['#', '$', '@'], skiprows=1, usecols=1)

	grids_size = np.loadtxt(script_dir + "/grid.dat", dtype=np.int64, comments=['#', '$', '@'], usecols=range(dim), max_rows=1)
	num_cell = np.prod(grids_size)

	ibz = np.loadtxt(script_dir + "/ibz.dat")

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

	alat=10.3310004

	for i in range(1,len(ibz)+1):
		i_path=f'./{path}/harmonic_{supercell}_dyn'
		q=qpoint(i_path,i)
		# q=[[0,0,0]]
		#print(i,len(q),q)
		print(np.dot(q,np.transpose(get_reduced_lattice())))
		o_path = f'./{path}/' + 'f2q_'
		output_new(q,i_path,o_path,i)
		print(f'file {o_path}{i} done')
		# exit()
