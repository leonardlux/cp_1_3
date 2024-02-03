
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
border_cells_per_side = 1   
# cells added at both spaces border to calculate at the borders
t_offset = 2               
# the first two point of times are needed to initalize the the simulation
# they are therefore not yet part of the distribution 
