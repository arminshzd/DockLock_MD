import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import numpy.ma as ma

# method to read free energy data created by WHAM
def read_F_data(fname):
    output = []
    inp_handle = open(fname, 'r')
    inp_content = inp_handle.readlines()
    inp_handle.close()
    inp_content = inp_content[1:]

    for i in inp_content:
        temp = []
        line = i.split()
        for j in line[:3]:
            temp.append(float(j))
        output.append(temp)

    return output

# Method to read trajectory files and super impose on the free energy surface
def read_traj_data(fname):
    r = []
    q = []
    inp_handle = open(fname, 'r')
    inp_content = inp_handle.readlines()
    inp_handle.close()

    for i in inp_content:
        line = i.split()
        r.append(float(line[1]))
        q.append(float(line[2]))

    return r, q


plt.figure(figsize=([20,10]), dpi=300)
data = read_F_data('1d_far.txt') # read free energy data

# Change this to read a trajectory file and superimpose on the free energy surface. Look at the loop below
traj_num = -1
if traj_num != -1:
    for i in range(traj_num):
        rtraj, qtraj = read_traj_data('./trajs/run_%g.traj' % i)
        if (rtraj[-1] < 2 and qtraj[-1] < 0):
	        plt.plot(rtraj, qtraj, 'ko-', markersize = 2)

# Plot the free energy surface
rsize, qsize = (90, 40)
r = np.linspace(1.05, 9.95, rsize) #80
q = np.linspace(-1.0725, 1.0725, qsize) #42
F = np.zeros((rsize,qsize))
for i in range(rsize):
    for j in range(qsize):
        F[i][j] = data[i*qsize+j][2]

F = F[:][:70] # Truncate the free energy surface to the region of interest
F = F.T
kbT = 1.380649e-23*300 # T = 300K
F = F*4180/6.022e23/kbT # convert to kT
[R, Q] = np.meshgrid(r[:70], q)
f_max = np.round(np.max(ma.masked_equal(F, np.inf)))
cmap = mpl.cm.Reds_r # change this to change the color map
norm = mpl.colors.Normalize(vmin=0, vmax=120)
plt.contourf(R, Q, F,  levels = 25, cmap = cmap, norm=norm)#'viridis' 'gist_ncar' levels=np.linspace(0, f_max)
cb = plt.colorbar()
plt.xlabel(r"$r_{COM}$ $(\AA)$", fontsize=25)
plt.ylabel(r"$q$", fontsize=25)
cb.set_label(r"Free Energy $(k_BT)$", fontsize=25)
plt.tight_layout()
plt.savefig("FES.png")

