import numpy as np

# Standard Parameters

# Parameters of our space:
x_min   = 0
x_max   = 1
dx      = 0.005
n_x     = int( ( x_max - x_min ) / dx + 1 )
# n_x is the amount of points for x on the spere (include 0 => +1)
# or [x_min,x_max] not (x_min,x_max]

# we just assume a square
y_min   = x_min
y_max   = x_max
dy      = dx
n_y     = int( ( y_max - y_min ) / dy + 1 )

assert dx == dy
# we used this assumption in the simulation

# time 
t_min   = 0
t_max   = 4
dt      = 0.002
n_t     = int (( t_max - t_min ) / dt + 1)

# Speed of wave:
c_c = 1 # c constant

# Task 3.4
sigma = 0.001
r_0_x = 0.5
r_0_y = 0.5

# Simulation paras:
border_cells_per_side = 2  
# cells added at both spaces border to calculate at the borders
t_offset = 2               
# the first two point of times are needed to initalize the the simulation
# they are therefore not yet part of the distribution 

class Config:
    def __init__(self) -> None:
        self.x_min  = x_min
        self.x_max  = x_max
        self.dx     = dx
        self.n_x    = n_x

        self.y_min   = y_min
        self.y_max   = y_max
        self.dy      = dy
        self.n_y     = n_y

        # time 
        self.t_min   = t_min
        self.t_max   = t_max
        self.dt      = dt
        self.n_t     = n_t

        # Speed of wave:
        self.c_c = c_c # c constant

        # Task 3.4
        self.sigma = sigma
        self.r_0_x = r_0_x
        self.r_0_y = r_0_y

        # Simulation paras:
        self.border_cells_per_side = border_cells_per_side  
        # cells added at both spaces border to calculate at the borders
        self.t_offset = t_offset 
        pass

    def recalc(self): 
        self.n_x     = int( ( self.x_max - self.x_min ) / self.dx + 1 )
        self.n_y     = int( ( self.y_max - self.y_min ) / self.dy + 1 )
        self.n_t     = int (( self.t_max - self.t_min ) / self.dt + 1)

        pass

    def test_stability(self):
        h = self.dx
        dt = self.dt
        c = self.c_c
        print("Numerical threshold anaylsis:")
        print(f"{h=}")
        print(f"{dt=}")
        print(f"{c=}")
        print("Equation: ")
        print(f"A = h/dt = {(h/dt):.3f}")
        print(f"B = c * √2 = {(c * np.sqrt(2)):.3f}")
        print(f"A > B: {h/dt > c * np.sqrt(2) }")
        pass