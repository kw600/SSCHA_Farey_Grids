import os,sys
import subprocess



output_path = 'target_666'


lte=['333','444']
mixing = True
#sscha final dyn matrix name
dyn_name=['dyn_pop2_','dyn_pop3_']

for i in range(len(lte)):
    if mixing and i==len(lte)-1:
        subprocess.run(["python", "mix.py", str(lte),str(dyn_name)])
    
    path  = lte[i]
    #generate the position file and force.dat for two smaller cells
    subprocess.run(["python", "generate_R.py", path])
    subprocess.run(["bash","q2r",path,dyn_name[i]])
    subprocess.run(["python", "normal_to_force.py", path])
    print('finish generating force.dat for '+path)


#genearate the position file  for the target large cell
subprocess.run(["python", "generate_R.py", output_path])
print('finish generating position file for '+output_path)

#read the input from two smaller cells and generate dynamical matrix for the target large cell to the output_path
subprocess.run(["python", "vcF2A.py", str(lte),output_path])
print('finish generating dynamical matrix for '+output_path)

# use fortran code to generate force.dat for the target large cell
os.chdir(output_path)
subprocess.run(["../a.out"])
os.chdir("..")
print('finish generating force.dat for '+output_path)

# transform the force.dat into qe format. Possible to use mpi4py to speed up the process for large supercell
subprocess.run(["python", "farey2qe.py", output_path])
print('finish generating qe format dynamical for '+output_path)
