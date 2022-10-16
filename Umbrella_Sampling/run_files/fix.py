import os
os.chdir("..")
kbT = 500*8.314*0.000239006 #kcal
rstart = 0.50
thstart = -1.1
rsteps = 71
thsteps = 45
rfconstant = 200#100
thfconstant = 100#50
whaminp_content = []
for i in range(thsteps):
    for j in range(rsteps):
        run_num = "{}-{}".format(i,j)
        centers = [rstart + j/20, thstart + i/20]
        fcs = [rfconstant, thfconstant]
        whaminp_content.append("out.{0}.traj\t{1}\t{2}\t{3}\t{4}\n".format(run_num, centers[0]*2, centers[1], fcs[0]*kbT/4, fcs[1]*kbT))

#for i in range(thsteps):
#    for j in range(rsteps):
#        run_num = "{}-{}".format(i,j)
#        os.chdir("./workingdir/{}".format(run_num))
#        print(run_num)
#        os.system("cp colvars.traj ../../output/out.{}.traj".format(run_num))
#        os.chdir("../..")

os.chdir("./output")
whaminp_handle = open("whaminput.pmf","w")
for i in whaminp_content:
    whaminp_handle.write(i)
whaminp_handle.close()
print("Done!")
