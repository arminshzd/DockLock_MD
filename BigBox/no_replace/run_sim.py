import os
import time


def run_MD(run_num):
        os.system(" cp ./dev_files/a.out ./dev_files/job.sh ./dev_files/init.coord ./runs/{}".format(run_num))
        os.chdir("./runs/{}".format(run_num))
        os.system("pwd")
        os.system("sbatch job.sh")
        os.chdir("../..")
        os.system("pwd")
        print("-"*150)

for i in range(100):
    run_num = "{}".format(i)
    os.system("mkdir ./runs/{}".format(run_num))
    run_MD(run_num)


print("Done!")
