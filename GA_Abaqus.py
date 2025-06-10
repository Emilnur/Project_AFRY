# Author: Balder Lohman Ranheim

import math
from deap import base, creator, tools, algorithms
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# Create a class for the Individual.
class MyIndividual:

    """ A class to represent an individual """

    def __init__(self, dynamic_attributes, static_attributes, geometry,stress_field):

        """
        Initialize an individual object.

        Parameters
        ----------
        dynamic_attributes (list): Attributes of the individual that can change.
        static_attributes (list): Attributes of the individual that are static.
        geometry (list): Geometry attributes of the individual.
        stress_field (list): The stress_field of the individual.

        """

        # Define the static attributes
        self.static_attributes = static_attributes

        # Define the dynamic attributes.
        self.dynamic_attributes = dynamic_attributes
        self.geometry = geometry
        self.stress_field = stress_field

        # Define the name for the individual.
        self.name = f"{dynamic_attributes[0],dynamic_attributes[1],dynamic_attributes[2]}"

    def calc_volume(self):
        
        """
        Calculates the volume of the individual based on the geometry and
        thicknesses.

        Return
        ----------
        vol_tot (float): The total volume of the individual.
        """
                
        # Calculates the total volume of material for the structure,
        # using the given geometry

        # Extract the geometry and thickness attributions.
        vR, nR, h, lN, rP = self.geometry
        thickness_vessel, thickness_Pad, thickness_Nozzle = self.dynamic_attributes

        # Calculate the outer radius.
        vR_out = vR + thickness_vessel
        nR_out = nR + thickness_Nozzle

        # Radiuses for the torus.
        R_rP = (rP + nR)/2
        r_rP = (rP - nR)/2 
         
        # Calculate the total volume of the vessel.
        vol_vessel = math.pi*h*(vR_out**2 - vR**2) - math.pi*thickness_vessel*nR**2 # Strip the material from where the nozzle is located.
        vol_nozzle = math.pi*lN*(nR_out**2 - nR**2)
        vol_pad = 2*math.pi**2**R_rP*r_rP**2 # Approximate as donout.

        vol_tot = vol_vessel + vol_nozzle + vol_pad

        return vol_tot,vol_vessel,vol_nozzle,vol_pad

    def calc_mass(self):

        """
        Calculates the mass of the individual based on the volume and density
        of the individual.

        Return
        ----------
        mass_tot (float): The total mass of the individual.
        """

        volume, _, _, _ = self.calc_volume()
        mass_tot = self.static_attributes['density'] * volume

        return mass_tot

    def compare_stresses(self):

        """
        Compares each of the stresses in the stress field for the given tolerance
        and return a value that either disqualifies the individual or return a neutral
        score.
        """

        for stress_field in self.stress_field:

            if max(stress_field) > self.static_attributes['stress_tol']:

                return np.inf

            else:

                pass

        return 0

    def fitness_func(self):

        """
        Calculate the fitness of the individual.
        
        Return
        ----------
        fitness (tuple): The fitness of the current individual.
        """        

        weight = self.calc_mass()
        stress_score = self.compare_stresses()

        mass_score = weight
        fitness = mass_score + stress_score # Might change based on other rewards.

        return fitness


def extract_stress_field(data_set):

    """
    Extracts the stress fields from each of the files in the data-set.

    Paramaters
    ----------
    data_set (list): Contains every path to the output files. 
    
    Return
    ----------
    selected_stresses (list): Contains every stress field from each of the
        files in the data-set.
    """

    selected_stresses = []
    for data_index,data in enumerate(data_set):

        df = pd.read_csv(data, delimiter=';')
        s11 = df['S11']
        s22 = df['S22']
        mises = df['VonMises']

        selected_stresses.append([s11, s22, mises])

    return selected_stresses


def extract_thickness_data(data_set):

    """
    Extracts each of the thicknesses from each of the files in the data-set.

    Paramaters
    ----------
    data_set (list): Contains every path to the output files. 
    
    Return
    ----------
    thicknesses_tot (list): Contains every defined thickness from the data-set.
    """

    thicknesses_tot = []
    for data in data_set:

        df = pd.read_csv(data)
        thickness_Vessel = df['Thickness_Vessel'] # Thickness of the vessel.
        thickness_Pad = df['Thickness_Pad'] # Thickness of the pad.
        thickness_Nozzle = df['Thickness_Nozzle'] # Thickness of the nozzle.

        for sim_index in range(len(thickness_Vessel)):

            thicknesses_tot.append([thickness_Vessel[sim_index]/1000,thickness_Pad[sim_index]/1000,thickness_Nozzle[sim_index]/1000])

    return thicknesses_tot

# Not finsihed function, fill in later.
def extract_geometry(data_set):

    """
    Extracts the geometry from each of the files in the data-set.

    Paramaters
    ----------
    data_set (list): Contains every path to the output files.
    
    Return
    ----------
    geometry (list): Contains the different geometry components.
    """

    geometry = []

    for data in data_set:

        df = pd.read_csv(data)

        vR = df['VR'] # Vessel radius.
        nR = df['NR'] # Nozzle radius.
        H = df['H'] # Height.
        lN = df['LN'] # Nozzle length.
        rP = df['RP'] # Reinforcement pad radius.

        for sim_index in range(len(vR)):

            geometry.append([vR[sim_index]/1000,nR[sim_index]/1000,H[sim_index]/1000,lN[sim_index]/1000,rP[sim_index]/1000])

    return geometry
    

def read_output_files(folder_path,result_file):

    """
    Read the output files defined for a specific folder.

    Paramaters
    ----------
    folder_path (string): A path to the folder for the output-files.
    
    Return
    ----------
    file_list (list): A list of the paths to the output files.
    """

    file_list = [str(file) for file in folder_path.rglob(f'*{result_file}')]

    return file_list

# Create the individual from the specifications in the initialization process.
def create_individual(individual_index,):

    """
    Creates an instance of an individual based on an index.

    Paramaters
    ----------
    individual_index (integer): An index for the current set of thickness
      (essentially what defines the individual). 
    
    Return
    ----------
    Individual (object): An instance of an individual.
    """

    dynamic_attributes = selected_thicknesses[individual_index]
    stress_field = selected_stress_fields[individual_index]
    geometry = selected_geometry[individual_index]
    
    return creator.Individual(
        dynamic_attributes=dynamic_attributes,static_attributes=static_attributes,geometry=geometry,stress_field=stress_field)

def evaluate(individual):

    """
    Evauates an individual by calling its corresponding fitness-value.

    Paramaters
    ----------
    indivudal (object): Instance of an individual.
    
    Return
    ----------
    fitness (tuple): The fitness-value for the individual.
    """

    return (individual.fitness_func(),) # Make sure it is a tuple

# Read all of the result files for the specified folder and collect them in a list.
path_directory = Path('GA_TestData')
data_set_results = read_output_files(path_directory,'pvessel*.csv')
data_set_geothick = read_output_files(path_directory,'doe*.csv')

# Extract necassary variables from the result files.
selected_geometry = extract_geometry(data_set_geothick)
selected_thicknesses = extract_thickness_data(data_set_geothick)
selected_stress_fields = extract_stress_field(data_set_results)

# Define the static input.
static_attributes = {

    'density': 7800, # [kg/m^3]
    'stress_tol': 423 # [Mpa]

}

## GA implementation.

# Create Type.
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", MyIndividual, fitness=creator.FitnessMin)

# Create the toolbox and register the components.
toolbox = base.Toolbox()
toolbox.register("evaluate", evaluate)
toolbox.register("select", tools.selTournament, tournsize=100) # The tournament size.

# Define the main program.
def main():

    # Create the population based on the amount of registered thicknesses (individuals).
    population = [create_individual(index) for index in range(len(selected_thicknesses))]

    # Register stats.
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # No mutation as for now.
    pop, logbook = algorithms.eaSimple(

        population, toolbox,
        cxpb=0.0, mutpb=0.0, ngen=10,
        stats=stats,
        halloffame=hof,
        verbose=True)

    return pop, logbook, hof

if __name__ == "__main__":
    
    pop, log, hof = main()
    print("Best individual is: %s\nwith fitness: %s" % (hof[0].name, hof[0].fitness))

    volume_tot,volume_vessel,volume_nozzle,volume_pad = hof[0].calc_volume()
    print(volume_vessel+volume_nozzle)
    print(volume_pad)
    print(hof[0].geometry)

    gen, avg, min_, max_ = log.select("gen", "avg", "min", "max")
    plt.plot(gen, avg, label="average")
    plt.plot(gen, min_, label="minimum")
    plt.plot(gen, max_, label="maximum")
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.legend(loc="lower right")
    plt.show()