import os
#import matplotlib
#matplotlib.use('GTK')
import matplotlib.pyplot as plt
import numpy as np

T = 500; #K
kb = 1.38065  #A^2.g/particle/K/s^2
bt = 1/T/kb
epsilon = 4e-8 #friction g/s
D = 2/bt/epsilon
smch = []
data_list = []
#os.chdir("./workingdir")
j = 14
k = 24
#os.chdir('./{}_{}'.format(j,k))
os.system('pwd')
traj_handle = open('rate_output.traj', 'r')
traj_content = traj_handle.readlines()
traj_handle.close()
for i in traj_content:
    temp_content = i.split()
    try:
        data_list.append((float(temp_content[0]), float(temp_content[1])))
    except:
        continue
docked = 0
locked = 0
incomplete = 0
for i in data_list: #(q, r)
    if i[1] < 2:
    	if i[0] > 0:
        	docked += 1
    	else:
        	locked += 1
    else:
        incomplete += 1

beta = (docked+locked)/len(data_list)
rate = 4*np.pi*D*j*beta/(1-(1-beta)*j/k)
rate *= 6.022e-4 
betainf = beta/(1-(1-beta)*j/k)
print("total number of runs: %d" % len(data_list))
print("total number of successful runs: %d (%%%f)" % (docked + locked ,beta*100))
print("rate constant of association: %f" % rate)
print("docked: %d\tlocked: %d" % (docked, locked))
print("percentage directly to locked: %f" %(locked/(docked+locked)*100))
per_docked = docked/(docked+locked)
rate_docked = per_docked*rate
print("Docking rate constant: %E" % rate_docked) # L/mol/s
