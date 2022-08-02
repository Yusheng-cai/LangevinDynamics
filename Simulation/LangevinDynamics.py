from openmm.app import *
import openmm.unit as unit
import openmm
import numpy as np
from tqdm import tqdm
from Simulation.Potential import BasePotential
import multiprocessing as mp

class SingleParticleSimulation:
    """
    Langevin Dynamics class that simulates the path of one particle under some potential 
    """
    def __init__(self, 
                    Potential: BasePotential,  
                    friction=10, 
                    temperature=300, 
                    timestep=2, 
                    mass=1.0, 
                    platform="CPU"):
                self.Potential_  = Potential
                self.friction_   = friction / unit.picosecond
                self.temperature_= temperature * unit.kelvin
                self.mass_       = mass * unit.dalton
                self.timestep_   = timestep * unit.femtosecond

                # Initialize integrator 
                self.integrator_ = openmm.LangevinIntegrator(temperature, friction, timestep)

                # Initialize the initial coordinates 
                self.coords_     = np.array([[0,0,0]])

                # add the particle to the potential
                self.Potential_.addParticle(0,[])

                # Add the system 
                self.system_     = openmm.System()
                self.system_.addParticle(self.mass_)
                self.system_.addForce(self.Potential_)

                # platform on which we are running the code 
                self.platform_   = openmm.Platform.getPlatformByName(platform)
                self.num_threads = str(mp.cpu_count())
                self.properties  = {"Threads" : self.num_threads}

                # set up context that stores the current complete state of a simulation
                self.context_    = openmm.Context(self.system_, self.integrator_, self.platform_, self.properties)
                self.context_.setPositions(self.coords_)
                self.context_.setVelocitiesToTemperature(self.temperature_)
    
    def __call__(self, iterations:int, print_frequency=50, outputfileName="output.dat"):
        #f = open(outputfileName, "w")
        #f.write("# x\ty\tz\tPE\tKE\tTE\tforce")
        for i in range(iterations):
            self.integrator_.step(1)

            # obtain the state object 
            state = self.context_.getState(getForces=True, getPositions=True, getVelocities=True, getEnergy=True)

            # obtain position
            position = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            print("Position = ", position)
            PE       = state.getPotentialEnergy()
            force    = state.getForces()
            print("Force from dw = " , force)
            KE       = state.getKineticEnergy()
            TE       = PE + KE

            # record in file
            #if (i+1) % print_frequency == 0:
            #    f.write("{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(position[0], position[1], position[2], PE, KE, TE))
                
        #f.close()

