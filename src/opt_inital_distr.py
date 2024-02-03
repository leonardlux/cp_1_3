import numpy as np
from . import config as c

def initial_sin(X,Y):
    return np.sin(np.pi*X)*np.sin(2*np.pi*Y)

def inital_exp(X,Y):
    return np.exp(-1/c.sigma * ( (X - c.r_0_x)**2 + (Y - c.r_0_y)**2))