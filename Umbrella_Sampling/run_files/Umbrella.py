import os
import time


def run_MD(centers, fcs, num_md_steps, run_num):
        rcenter = centers[0]
        thcenter = centers[1]
        rfconstant = fcs[0]
        thfconstant = fcs[1]
        print('rcenter: {}\nthcenter: {}'.format(rcenter,thcenter))
        col_handle = open("main.c","r")
        col_content = col_handle.readlines()
        col_handle.close()
        rcenter_ind = col_content.index("double r0U = 0;\n")
        thcenter_ind = col_content.index("double q0U = 0;\n")
        num_steps_ind = col_content.index("unsigned long long num_runs = 0;\n")
        rfc_ind = rcenter_ind-1
        thfc_ind = thcenter_ind-1
        new_col_content = col_content
        new_col_content[rcenter_ind] = "double r0U = {};\n".format(rcenter)
        new_col_content[thcenter_ind] = "double q0U = {};\n".format(thcenter)
        new_col_content[rfc_ind] = "double krU = {};\n".format(rfconstant)
        new_col_content[thfc_ind] = "double kqU = {};\n".format(thfconstant)
        new_col_content[num_steps_ind] = "unsigned long long num_runs = {};\n".format(num_md_steps)
        col_temp_handle = open("main_temp.c","w")
        for j in new_col_content:
                col_temp_handle.write(j)
        col_temp_handle.close()
        os.system(" cp main_temp.c job.sh init.coord ../workingdir/{}".format(run_num))
        os.chdir("../workingdir/{}".format(run_num))
        os.system("pwd")
        os.system("gcc main_temp.c -lm -lgsl -lgslcblas -std=gnu99; sbatch job.sh")
        os.chdir("../../run_files")
        os.system("pwd")
        print("-"*150)

rstart = 0.50
thstart = -1.1
rstepss = 0
rstepse = 71
thstepss = 0
thstepse = 45
rfconstant = 200#100
thfconstant = 100#50
num_md_steps  = 1000
single_runtime = 120
active_runs = 0
for i in range(thstepss, thstepse):
    for j in range(rstepss, rstepse):
        run_num = "{}-{}".format(i,j)
        os.system("mkdir ../workingdir/{}".format(run_num))
        centers = [rstart + j/20, thstart + i/20]
        fcs = [rfconstant, thfconstant]
        run_MD(centers, fcs, num_md_steps, run_num)


print("Done!")
