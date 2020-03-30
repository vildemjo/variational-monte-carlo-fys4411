# Common imports
import os
# Where to save the figures and data files

from numpy import log2, zeros, mean, var, sum, loadtxt, arange, array, cumsum, dot, transpose, diagonal, sqrt
from numpy.linalg import inv

import pandas as pd
from pandas import DataFrame

def block(x):
    # preliminaries
    n = len(x)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)

    # estimate the auto-covariance and variances 
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
   
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")

    return mu, s[k]/2**(d-k)


def calculate_mean_squared(x):
    mean_squared = 0.0
    for i in range(len(x)):
        mean_squared += x[i]*x[i]
    mean_squared /= float(len(x))
    return mean_squared

def calculate_uncor_var(x):
    uncor_var = 0.0
    x0 = mean(x)
    for i in range(len(x)):
        uncor_var += (x[i] - x0)**2
    uncor_var /= float(len(x))
    return uncor_var

def exact_energy(alpha, N):
    return (0.5*alpha + 1/(8.0*alpha))*N*3

# exercise b

# DATA_ID = "../Output//exercise_e//allEnergies"

# Ns = ["50"]#, "50", "100", "500"]#, "100"]#, "500"]

# As = ["50"]#[ "35", "40", "45", "55", "60", "65"]

# for j in range(len(Ns)):
#     print "N: ", Ns[j]
#     print "$alpha$: & $left< E_L right>$:& $E_{exact}$ & $sigma_B$ & $sigma$"+ "\\" + "\\"
    
#     for a in range(len(As)):    
        
#         def data_path(dat_id):
#             return os.path.join(DATA_ID, dat_id)

#         infile = open(data_path("analytical_3d_%sp_alpha_%s_MC_20_energy.txt"%(Ns[j],As[a])),'r')

#         x = loadtxt(infile)
#         x = x[:int(2**19)]

#         (mu, variance) = block(x) 
#         std = sqrt(variance)

#         uncor_std = sqrt(calculate_uncor_var(x)) # sqrt(calculate_mean_squared(x) - mu*mu)
        
#         alpha = float("0.%s"%As[a])
#         N = float(Ns[j])

#         print "0.%s & %.5f & %.5f & %.5f & %.5f & %.5f"%(As[a], mu, exact_energy(alpha,N), abs(exact_energy(alpha,N)-mu) , std, uncor_std)+ "\\"+ "\\"



# -------------------------------------------

DATA_ID = "../Output//"

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

infile = open(data_path("exercise_f/50p_3d_post_gradient_descent_interacting_energy.txt"),'r')

x = loadtxt(infile)

(mu, variance) = block(x) 
std = sqrt(variance)

uncor_std = sqrt(calculate_uncor_var(x)) # sqrt(calculate_mean_squared(x) - mu*mu)

print " 0.483195 & %.5f & %.5f & %.5f"%( mu, std, uncor_std)+ "\\"+ "\\"



