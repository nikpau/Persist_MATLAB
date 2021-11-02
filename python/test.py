import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import numpy as np


fig = plt.figure()

gs = fig.add_gridspec(3,3)

f_ax1 = fig.add_subplot(gs[:, :-1])
f_ax1.set_title('gs[0, :]')
f_ax2 = fig.add_subplot(gs[0, -1])
f_ax2.set_title('gs[1, :-1]')
f_ax3 = fig.add_subplot(gs[1, -1])
f_ax3.set_title('gs[1:, -1]')
f_ax4 = fig.add_subplot(gs[2, -1])
f_ax4.set_title('gs[-1, 0]')

import time


a=1
b=4
c=61
d=5
e=8
f=9

def rand_calc_global(a,b,c,d,e,f):

    x = a + b**(np.sqrt(c/f))
    y = a*b*c//d **(e*f)

    z = x**(np.sqrt(y))

    return z

def rand_calc_local():

    a=1
    b=4
    c=61
    d=5
    e=8
    f=9

    x = a + b**(np.sqrt(c/f))
    y = a*b*c//d **(e*f)

    z = x**(np.sqrt(y))

    return z

# Test loop

t_start = time.time()

def main(mode):


    for _ in range(100000):
        if mode == "global":    
            p = rand_calc_global(a,b,c,d,e,f)

            t_end = time.time()

            print(f"Time for", mode, ":", t_end-t_start )

        else:
            p = rand_calc_local()

            t_end = time.time()
            
            print(f"Time for", mode, ":", t_end-t_start )

main("global")