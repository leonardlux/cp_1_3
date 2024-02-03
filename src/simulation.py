import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ipywidgets import interact

# custom files
from . import config as c
from .opt_inital_distr import initial_sin

class Simulation:
    def __init__(self,func_u_0):
        # border cells for the edge of the calculation
        self.bc = c.border_cells_per_side + 1
        # +1 to use the -1 for the last relevant index
        self.x_dim   = c.n_x + (self.bc * 2)
        self.y_dim   = c.n_y + (self.bc * 2)
        self.t_dim   = c.n_t + c.t_offset

        # initalize U 
        # array of matricies
        self.U = np.zeros(( self.t_dim , self.x_dim , self.y_dim ))

        # arrays
        x              = np.linspace( c.x_min, c.x_max, c.n_x ) 
        y              = np.linspace( c.y_min, c.y_max, c.n_y ) 
         # meshgrid for broadcasting calculation
        self.X, self.Y = np.meshgrid( x, y)                    

        # t is just needed for plotting
        self.t   = np.linspace(
                (c.t_min - c.t_offset * c.dt),
                c.t_max,
                c.n_t + c.t_offset,
            ) 
        # adjusting, so that the simulation starts at t_min
        
        # ignore the bordercells, they will be set to 0
        # this is already an implicit boundary condition
        self.U[ 0, self.bc:-self.bc, self.bc:-self.bc ] = func_u_0( self.X , self.Y )

        # calculate the constant 
        self.alpha = c.dt**2 * c.c_c**2 / c.dx**2
        # here we use the assumption that c.dx == c.dy
        
        # Set up t=1 second condition (du/dt = 0 for t = 0 => t=1 calculated differently)
        self.U[ 1, self.bc:-self.bc, self.bc:-self.bc ] = \
            + self.U[ 0, self.bc:-self.bc, self.bc:-self.bc] \
            + 1/2 * self.alpha * (   
                    (                 # d^2/dx^2
                        +     self.U[ 0, self.bc +1 :-self.bc +1, self.bc    :-self.bc   , ] 
                        +     self.U[ 0, self.bc -1 :-self.bc -1, self.bc    :-self.bc   , ]
                    ) + (                            # d^2/dy^2
                        +     self.U[ 0, self.bc    :-self.bc   , self.bc +1 :-self.bc +1, ] 
                        +     self.U[ 0, self.bc    :-self.bc   , self.bc -1 :-self.bc -1, ]
                    )  # part of both:
                    - 4 *     self.U[ 0, self.bc    :-self.bc   , self.bc    :-self.bc   , ] 
                )
        
    def run(self):
        # iterate through time
        # only start after the inital conditions (c.t_offset)
        for i_t in range( c.t_offset, self.t_dim ):
            self.U[ i_t, self.bc:-self.bc, self.bc:-self.bc ] = \
                + 2 * self.U[ i_t-1, self.bc:-self.bc, self.bc:-self.bc ] \
                - 1 * self.U[ i_t-2, self.bc:-self.bc, self.bc:-self.bc ] \
                + self.alpha * (   
                    (                 # d^2/dx^2
                        +     self.U[ i_t-1, self.bc +1 :-self.bc +1, self.bc    :-self.bc   , ] 
                        +     self.U[ i_t-1, self.bc -1 :-self.bc -1, self.bc    :-self.bc   , ]
                    ) + (                                      # d^2/dy^2
                        +     self.U[ i_t-1, self.bc    :-self.bc   , self.bc +1 :-self.bc +1, ] 
                        +     self.U[ i_t-1, self.bc    :-self.bc   , self.bc -1 :-self.bc -1, ]
                    )  # part of both x and y:
                    - 4 *     self.U[ i_t-1, self.bc    :-self.bc   , self.bc    :-self.bc   , ] 
                )
        #TODO: do i need to implement the border condition here? Is it already part of it ?
    
    def plot(self,z_lim=1.1):
        fig  =  plt.figure(figsize = (12, 5))
        ax0  = fig.add_subplot(131, projection = '3d')
        ax1  = fig.add_subplot(132, projection = '3d')
        ax2  = fig.add_subplot(133, projection = '3d')

        axs = [ ax0, ax1, ax2 ]
        times = [
            0,              # start point
            int(c.n_t/2),   # half way through
            c.n_t -1,          # end of simulation
        ]
        # adjust for offset
        if c.t_offset !=0 :
            times = [t_i + c.t_offset for t_i in times]

        # plot everything
        for ax, t_i in zip(axs,times): 
            ax.plot_surface(self.X, self.Y, self.U[ t_i, self.bc:-self.bc, self.bc:-self.bc,])
            ax.set_zlim(-z_lim, z_lim)


        plt.tight_layout()
        plt.show()
    
    def plot_interactive(self,t_i_steps=100,z_lim=1.1):
        def plotter(t):
            t = int(t)
            fig  =  plt.figure(figsize = (12, 5))
            ax  = fig.add_subplot(projection = '3d')
            ax.plot_surface(self.X, self.Y, self.U[ t, self.bc:-self.bc, self.bc:-self.bc,])
            ax.set_zlim(-z_lim, z_lim)

            plt.title(f"at $t={t*c.dt}$")
        interact(plotter, t = (0, c.n_t + c.t_offset -1, t_i_steps));


if __name__ == "__main__":
    sim = Simulation(initial_sin)
    sim.run()
    sim.plot()