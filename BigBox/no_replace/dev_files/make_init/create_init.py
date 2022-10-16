import numpy as np

# Configuration:
no_mono = 50 #Number of monomers
fib_len = 10  #Fibril length
box_size = 60
pad = 5 #Padding around particles


# Initiations:
pos = np.zeros((2*(no_mono+fib_len), 3))
fib_cntr = np.zeros((fib_len, 3)) #center of fibrils in monomer

# Putting fibril in:

if fib_len%2==0:
    y_cntr = np.arange(-1*fib_len/2, fib_len/2)*0.707
else:
    y_cntr = np.arange(-1*(fib_len-1)/2, (fib_len-1)/2+1)*0.707
for i in range(fib_len):
    if i%2==0:
        if (2*i//4)%2==0:
            pos[2*i] = np.array([0, y_cntr[i], -0.5])
            pos[2*i +1] = np.array([0, y_cntr[i], 0.5])
        else:
            pos[2*i] = np.array([0, y_cntr[i], 0.5])
            pos[2*i +1] = np.array([0, y_cntr[i], -0.5])
    else:
        if (2*i//4)%2==0:
            pos[2*i] = np.array([0.5, y_cntr[i], 0])
            pos[2*i +1] = np.array([-0.5, y_cntr[i], 0])
        else:
            pos[2*i] = np.array([-0.5, y_cntr[i], 0])
            pos[2*i +1] = np.array([0.5, y_cntr[i], 0])
    fib_cntr[i] = np.array([0.5, y_cntr[i], 0])

# putting monomer in the box

flag1 = True
flag2 = True
flag3 = True
cnt = 0
length = 0
centers = np.zeros((no_mono,3))
while(flag1):
    rndcent = np.random.choice([1,-1], 3)*np.random.rand(3)*box_size
    for i in range(fib_len):
        if np.linalg.norm(rndcent-pos[i]) <= 0.5+pad:
            flag2 = False
            break
    if flag2:
        for i in range(fib_len+1, fib_len+1+length):
            if np.linalg.norm(rndcent-pos[i]) <= 1:
                flag3 = False
                break
    if flag3 or flag2:
        centers[length] = rndcent
        length += 1
    if length == no_mono:
        flag1 = False

for i in range(no_mono):
    sinth = np.random.choice([1,-1]) * np.random.rand()
    costh = np.random.choice([1,-1]) * np.sqrt(1-sinth**2)
    sinphi = np.random.choice([1,-1]) * np.random.rand()
    cosphi = np.random.choice([1,-1]) * np.sqrt(1-sinphi**2)
    rvec = np.array([sinphi*costh/2, sinphi*sinth/2, cosphi/2])
    pos[2*(i+fib_len)] = centers[i]+rvec
    pos[2*(i+fib_len)+1] = centers[i]-rvec


# Writing to file
fhandle = open("init.coord", 'w')
fhandle.write(str((no_mono+fib_len)))
fhandle.write("\n")
for i in range(2*(no_mono+fib_len)):
    fhandle.write(str(pos[i][0]))
    fhandle.write("\n")
for i in range(len(pos)):
    fhandle.write(str(pos[i][1]))
    fhandle.write("\n")
for i in range(len(pos)):
    fhandle.write(str(pos[i][2]))
    fhandle.write("\n")

fhandle.close()
