import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ipywidgets import interact

# custom files
from .opt_inital_distr import initial_sin

class Simulation:
    def __init__(self,func_u_0,config):
        self.c = config
        # border cells for the edge of the calculation
        self.bc = self.c.border_cells_per_side
        # +1 to use the -1 for the last relevant index
        self.x_dim   = self.c.n_x + (self.bc * 2)
        self.y_dim   = self.c.n_y + (self.bc * 2)
        self.t_dim   = self.c.n_t + self.c.t_offset

        # initalize U 
        # array of matricies
        self.U = np.zeros(( self.t_dim , self.x_dim , self.y_dim ))

        # arrays
        x              = np.linspace( self.c.x_min, self.c.x_max, self.c.n_x ) 
        y              = np.linspace( self.c.y_min, self.c.y_max, self.c.n_y ) 
         # meshgrid for broadcasting calculation
        self.X, self.Y = np.meshgrid( x, y)                    

        # t is just needed for plotting
        self.t   = np.linspace(
                (self.c.t_min - self.c.t_offset * self.c.dt),
                self.c.t_max,
                self.c.n_t + self.c.t_offset,
            ) 
        # adjusting, so that the simulation starts at t_min
        
        # ignore the bordercells, they will be set to 0
        # this is already an implicit boundary condition
        self.U[ 0, self.bc:-self.bc, self.bc:-self.bc ] = func_u_0( self.X , self.Y )

        # calculate the constant 
        self.alpha = self.c.dt**2 * self.c.c_c**2 / self.c.dx**2
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
        for i_t in range( self.c.t_offset, self.t_dim ):
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
        # Due to the limited scope of each calculation in the x, and y space (bc:-bc)
        # we we already fullfilled the boundary condtions.
    
    def plot(self,z_lim=1.1,analy=False, times = [],save_to = ""):
        def analy_sol(x,y,t):
            return np.cos(np.sqrt(5)*np.pi*self.c.dt* t) * np.sin(np.pi * x) * np.sin(2* np.pi * y)

        fig  =  plt.figure(figsize = (12, 5))
        ax0  = fig.add_subplot(131, projection = '3d')
        ax1  = fig.add_subplot(132, projection = '3d')
        ax2  = fig.add_subplot(133, projection = '3d')

        axs = [ ax0, ax1, ax2 ]
        if times == []:
            times = [
                0,              # start point
                int(self.c.n_t/2),   # half way through
                self.c.n_t -1,          # end of simulation
            ]
        # adjust for offset
        if self.c.t_offset !=0 :
            times = [t_i + self.c.t_offset for t_i in times]

        # plot everything
        for ax, t_i in zip(axs,times):
            ax.set_title(f"t = {self.c.dt*(t_i - self.c.t_offset)}")
            ax.plot_surface(
                self.X, 
                self.Y, 
                self.U[ t_i, self.bc:-self.bc, self.bc:-self.bc,],
                label="numerical",
                )
            ax.set_zlim(-z_lim, z_lim)
            if analy: 
                ax.plot_surface(
                    self.X, 
                    self.Y, 
                    analy_sol(self.X,self.Y,t_i),
                    label="analytical",
                    alpha= 0.7,
                    )
                plt.legend()

            ax.set_xlabel("$x$")
            ax.set_ylabel("$y$")
            ax.set_zlabel("$u(x,y,t)$")

        if save_to =="":
            plt.show()
        else:
            plt.savefig(save_to)
    
    def plot_interactive(self,t_i_steps=100,z_lim=1.1,analy=False,diff=False):
        def analy_sol(x,y,t):
            return np.cos(np.sqrt(5)*np.pi*self.c.dt* t) * np.sin(np.pi * x) * np.sin(2* np.pi * y)
        def plotter(t_i):
            t_i = int(t_i)
            fig  =  plt.figure(figsize = (12, 5))
            ax  = fig.add_subplot(1,2,1,projection = '3d')
            
            ax.plot_surface(self.X, self.Y, self.U[ t_i, self.bc:-self.bc, self.bc:-self.bc,])
            ax.set_zlim(-z_lim, z_lim)
            if analy: 
                ax.plot_surface(
                    self.X, 
                    self.Y, 
                    analy_sol(self.X,self.Y,t_i),
                    label="analytical",
                    alpha= 0.7,
                    )
            if diff:

                ax2  = fig.add_subplot(1,2,2,projection = '3d')
                ax2.plot_surface(
                    self.X, 
                    self.Y, 
                    (self.U[ t_i, self.bc:-self.bc, self.bc:-self.bc,] - analy_sol(self.X,self.Y,t_i))/np.max(self.U[ t_i, self.bc:-self.bc, self.bc:-self.bc,])*100 ,
                    )
                
            

            plt.title(f"at $t={t_i*self.c.dt}$")
        interact(plotter, t_i = (0, self.c.n_t + self.c.t_offset -1, t_i_steps));


if __name__ == "__main__":
    sim = Simulation(initial_sin)
    sim.run()
    sim.plot()