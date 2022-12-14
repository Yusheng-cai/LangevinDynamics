from openmm import unit
from openmm import openmm
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
                self.integrator_ = openmm.LangevinIntegrator(self.temperature_, self.friction_, self.timestep_)

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
        f = open(outputfileName, "w")
        f.write("# x(nm)\ty(nm)\tz(nm)\tPE(kJ/mol)\tKE(kJ/mol)\tTE(kJ/mol)\n")
        for i in tqdm(range(iterations)):
            self.integrator_.step(1)

            # obtain the state object 
            state = self.context_.getState(getForces=True, getPositions=True, getVelocities=True, getEnergy=True)

            # obtain position
            position = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
            PE       = state.getPotentialEnergy() / unit.kilojoule_per_mole
            force    = state.getForces() 
            KE       = state.getKineticEnergy() /unit.kilojoule_per_mole
            TE       = PE + KE

            # record in file
            if (i+1) % print_frequency == 0:
                f.write("{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(position[0,0], position[0,1], position[0,2], PE, KE, TE))
                
        f.close()

