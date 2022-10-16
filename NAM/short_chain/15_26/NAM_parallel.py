#!/share/apps/anaconda/anaconda3/bin/python3
from ctypes import CDLL
from joblib import Parallel, delayed

# Inputs
so_file = "/home/armin/vwall/NAM_pl/short_chain/15_26/main.so"
num_steps = int(1e10)
num_runs = 20000
num_cores = 96

# Function
def run_MD(num_steps, run_num):
    main_func = CDLL(so_file)
    md_return = main_func.run(num_steps, run_num)
    return md_return

# Parallelization

output = Parallel(n_jobs=num_cores)(delayed(run_MD)(num_steps, n) for n in range(num_runs))
