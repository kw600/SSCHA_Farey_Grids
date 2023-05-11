import os,sys
import subprocess



output_path = ['target_121212']

try:
	os.mkdir(output_path)
except:
	pass

lte=['333','444']
#sscha final dyn matrix name
dyn_name=['dyn_pop3_','dyn_pop3_']

for i in range(len(lte)):
	path  = lte[i]
	#generate the position file and force.dat for two smaller cells
	subprocess.run(["python", "generate_R.py", path])
	subprocess.run(["q2r",path,dyn_name[i]])
	subprocess.run(["python", "normal_to_force.py", path])
	print('finish generating force.dat for '+path)

#genearate the position file  for the target large cell
subprocess.run(["python", "generate_R.py", output_path])
print('finish generating position file for '+output_path)

#read the input from two smaller cells and generate dynamical matrix for the target large cell to the output_path
subprocess.run(["python", "vcF2A.py", path,output_path])
print('finish generating dynamical matrix for '+output_path)

# use fortran code to generate force.dat for the target large cell
os.chdir(output_path)
subprocess.run(["./a.out"])
os.chdir("..")
print('finish generating force.dat for '+output_path)

# transform the force.dat into qe format. Possible to use mpi4py to speed up the process for large supercell
subprocess.run(["python", "farey2qe.py", path,output_path])
print('finish generating qe format dynamical for '+output_path)
