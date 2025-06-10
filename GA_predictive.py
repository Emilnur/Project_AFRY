# Author: Balder Lohman Ranheim

import math
import random
from deap import base, creator, tools
from itertools import permutations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from Prediction_forGA import stress_eval
import time

# Time for the prediction model.
time_instances_tot = []

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
        vR_in = vR
        vR_out = vR_in + thickness_vessel
        nR_out = nR
        nR_in = nR_out - thickness_Nozzle

        # Radiuses for the pad.
        r_in_pad = nR - thickness_Nozzle/2 
        R_out_pad = rP

        # Calculate the total volume of the vessel.
        vol_vessel = math.pi*h*(vR_out**2 - vR_in**2) - math.pi*thickness_vessel*nR_in**2 # Strip the material from where the nozzle is located.
        vol_nozzle = math.pi*lN*(nR_out**2 - nR_in**2)
        vol_pad = thickness_Pad*(math.pi*(R_out_pad**2 - r_in_pad**2))

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

        volume_tot,_,_,_ = self.calc_volume()
        mass_tot = self.static_attributes['density'] * volume_tot

        return mass_tot

    def compare_stresses(self):

        """
        Compares each of the stresses in the stress field for the given tolerance
        and return a value that either disqualifies the individual or return a neutral
        score.
        """
        
        for max_stress in self.stress_field:

            if max_stress > self.static_attributes['stress_tol']:

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

# =======================================================================================

def predictive_model(thicknesses):

    """
    Provides the predicted max stresses for a given set of thicknesses.
    
    Paramaters
    ----------
    thicknesses (list): A set of thicknesses corresponding to an individual.
    
    Return
    ----------
    stress_max (list): The max values of the stress fields.
    """     

    # Extract each thickness and feed them into the predictive model.
    print(thicknesses)
    t1 = thicknesses[0]*1e3
    t2 = thicknesses[1]*1e3
    t3 = thicknesses[2]*1e3
    time1 = time.time()
    stress_max = stress_eval(t1,t2,t3)
    time2 = time.time()
    print(stress_max[0])
    time_tot = time2-time1
    print("The prediction process took: %s s" %time_tot)
    time_instances_tot.append(time_tot)

    return stress_max

# ----------------------------------------------------------------------------------------------

def thickness_gen_init(dyn_range, set_size, n_sets, initial_decimals):
    """
    Generate a specified number of unique permutations of thickness values
    within a given interval. Starts with a given decimal precision and 
    automatically increases precision and granularity as needed.

    Parameters
    ----------
    dyn_range (tuple): The range for the variation of the dynamic variables.
    set_size (int): Number of thicknesses in each set.
    n_sets (int): Number of unique sets to generate.
    initial_decimals (int): Starting number of decimal places.

    Return
    ----------
    thickness_sets (list): Each set representing a unique thickness set.
    """
    start, end = dyn_range

    decimals = initial_decimals
    n_points = set_size
    while True:
        
        while True:
            
            # 1. Generate values
            raw_values = np.linspace(start, end, n_points)

            # 2. Round to current precision
            values = np.round(raw_values, decimals)

            # 3. Keep only unique values that fall within the original interval
            values = np.unique(values)
            values = values[(values >= start) & (values <= end)]

            # Skip if not enough distinct values
            if len(values) < set_size:
                n_points += 1
                continue

            # 4. Generate permutations
            all_perms = list(permutations(values, set_size))

            if len(all_perms) >= n_sets:
                thickness_sets = [list(p) for p in random.sample(all_perms, n_sets)]
                return thickness_sets
            
            n_points += 1

            if n_points > 500:
                break  # fallback to increasing precision
        print("Decimals has incresed")
        decimals += 1
        n_points = set_size  # reset granularity

# ----------------------------------------------------------------------------------------------

def thickness_gen(num_dynamic,dyn_range):

    """
    Defines a set of thicknesses that all are within the defined range.
    
    Paramaters
    ----------
    num_dynamic (integer): The number of dynamic variables for each individuals.
    dyn_range (tuple): The range for the variation of the dynamic variables.

    Return
    ----------
    thicknesses (list): The defined thicknesses.
    """     

    # Extract the minimum and maximum value from the given range and generate
    # a thickness set that is within the upper and lower limit.
    min_dynamic, max_dynamic = dyn_range
    thicknesses = [random.uniform(min_dynamic,max_dynamic) for _ in range(num_dynamic)]

    return thicknesses

# ----------------------------------------------------------------------------------------------

def init_process(num_individuals,num_dynamic,dyn_range):

    """
    Initializes the process by creating a specified number of thickness sets
    corresponding to their separate max stresses.

    Paramaters
    ----------
    num_individuals (integer): A specification of the number of individuals.
    num_dynamic (integer): The number of dynamic variables for each individuals.
    dyn_range (tuple): The range for the variation of the dynamic variables.

    Return
    ----------
    thicknesses_tot (list): All of the thickness sets for the defined number of individuals.
    max_stresses_curr (list): All of the max stresses for the defined number of individuals.
    """

    max_stresses_tot = []

    # Generate all of the thickness sets.
    thickness_sets_tot = thickness_gen_init(dyn_range, num_dynamic, num_individuals, initial_decimals=6)
    for individual_ind in range(num_individuals):

        # Extract the max stresses.
        max_stresses_curr = predictive_model(thickness_sets_tot[individual_ind])

        # Store the current dynamic values and the max stresses for that correspond to the
        # current individual.
        max_stresses_tot.append(max_stresses_curr)

    return [thickness_sets_tot,max_stresses_tot]

# ----------------------------------------------------------------------------------------------

# Create the individual from the specifications in the initialization process.
def create_individual(individual_index,thicknesses,stress_field):

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

    dynamic_attributes = thicknesses[individual_index]
    stress_field = stress_field[individual_index]
    
    return creator.Individual(
        dynamic_attributes=dynamic_attributes,static_attributes=static_attributes,geometry=geometry,stress_field=stress_field)

# ----------------------------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------------------------

def mutate(individual,static_attr,geom,dyn_range):

    """
    Mutates an individual by changing one of the genes (one of the dynamic variables).

    Paramaters
    ----------
    indivudal (object): Instance of an individual.
    static_attr (dictionary): A definiton of the static attributes.
    geom (list): The attributes for the geometry.
    dyn_range (tuple): The range for the variation of the dynamic variables.
    
    Return
    ----------
    mutated_individual (object): Instance of a mutated individual.
    """

    # Extract the value corresponding to the thickness that is to be mutated.
    thickness_set = individual.dynamic_attributes
    thickness_index = random.randrange(len(thickness_set))
    thickness_value = thickness_set[thickness_index]

    # Extract the minimum and the maximum value in the dynamic range.
    min_thickness,max_thickness = dyn_range

    # Mutate the selcted thickness.
    mutated_thickness = thickness_value
    while mutated_thickness == thickness_value:

        mutated_thickness = random.uniform(min_thickness,max_thickness)

    # Update the thickness-set with the mutated value.
    thickness_set[thickness_index] = mutated_thickness

    # Calculate the max-stresses corresponding to the new thickness-set.
    max_stresses_mutated = predictive_model(thickness_set)

    # Create an instance of the mutated individual.
    mutated_individual = MyIndividual(dynamic_attributes=thickness_set,static_attributes=static_attr,geometry=geom,stress_field=max_stresses_mutated)

    return mutated_individual,

# ----------------------------------------------------------------------------------------------

def crossover(parent_1,parent_2,static_attr,geom,dyn_range):

    """
    Creates 2 child from 2 parents by mixing the parents genes (dynamic attributes).

    Paramaters
    ----------
    parent_1 (object): Instance of an individual.
    parent_2 (object): Instance of an individual.
    static_attr (dictionary): A definiton of the static attributes.
    geom (list): The attributes for the geometry.
    dyn_range (tuple): The range for the variation of the dynamic variables.
    
    Return
    ----------
    child1 (object): Instance of an individual.
    child2 (object): Instance of an individual.
    """

    # Extract the minimum and maximum value from the dynamic range.
    min_dyn,max_dyn = dyn_range

    # Create a child instance for every parent.
    child_genes_tot = []
    for _ in range(2):

        # Crossover the genes and keep the new genes within bounds.
        child_genes = []
        for par_ind in range(len(parent_1.dynamic_attributes)):
            
            # Scale factor for the genes, so child1 =! child2.
            sf = random.random()

            child_gene = min(max_dyn,max(min_dyn,np.mean([parent_1.dynamic_attributes[par_ind],parent_2.dynamic_attributes[par_ind]])*sf))
            child_genes.append(child_gene)
        
        child_genes_tot.append(child_genes)

    # Calculate the max-stresses corresponding to the new thickness-set.
    child1_max_stresses = predictive_model(child_genes_tot[0])
    child2_max_stresses = predictive_model(child_genes_tot[1])

    # Create an instance of the child.
    child1 = MyIndividual(dynamic_attributes=child_genes,static_attributes=static_attr,geometry=geom,stress_field=child1_max_stresses)
    child2 = MyIndividual(dynamic_attributes=child_genes,static_attributes=static_attr,geometry=geom,stress_field=child2_max_stresses)

    return child1,child2

# ========================================================================================================================================

# Define the size of the population, number of dynamic attributes and the range of the dynamic attributes.
number_of_individuals = 20000
number_of_dynamic_variables = 3
dynamic_range = (0.009,0.015)

# Create data with given definitions.
thicknesses_init,max_stresses_init = init_process(number_of_individuals,number_of_dynamic_variables,dynamic_range)

print("Init-process initalized")
# Define the static input.
static_attributes = {

    'density': 7800, # [kg/m^3]
    'stress_tol': 338.4 # [Mpa], with the capture of 20% error.

}

# The geometry for the given instance.
geometry = [2,0.508,2.6,0.3,0.66]

# Predefine list for tracking best individual.
best_individuals = []

# ========================================================================================================================================

# Create Type.
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", MyIndividual, fitness=creator.FitnessMin)

# Create the toolbox and register the components.
toolbox = base.Toolbox()
toolbox.register("evaluate", evaluate)
toolbox.register("select", tools.selTournament, tournsize=100) # The tournament size.
toolbox.register("mate", crossover, static_attr=static_attributes,geom=geometry,dyn_range=dynamic_range)
toolbox.register("mutate", mutate, static_attr=static_attributes,geom=geometry,dyn_range=dynamic_range)

# The genetic algorithm.
def main():

    """
    The genetic algorithm, modified to fit the continous updates of the stress fields
    for new dynamic attributes through mutation or crossover.
    
    Return
    ----------
    population (list): The population made up of individuals.
    logbook (object): Stores log entries.
    hof (object): The best individuals from all generations.
    """

    # Create the population based on the amount of registered thicknesses (individuals).
    population = [create_individual(index,thicknesses_init,max_stresses_init) for index in range(len(thicknesses_init))]    

    # Register stats.
    hof = tools.HallOfFame(1)
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)

    # Define the number of generations, crossover probability and the mutation probability.
    ngen = 5
    cxpb = 0.01
    mutpb = 0.01

    # Track the progress by registrering stats.
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + stats.fields

    # Initial evaluation.
    fitnesses = map(toolbox.evaluate, population)
    print(fitnesses)
    for ind, fit in zip(population, fitnesses):
        ind.fitness.values = fit

    hof.update(population)
    record = stats.compile(population)
    logbook.record(gen=0, nevals=len(population), **record)
    print(logbook.stream)

    # Save the best indivudal for the current generation.
    best_individuals.append(hof[0].fitness)
    
    # Go through each generation.
    for gen in range(1, ngen+1):

        offspring = []
        # Check if population has been replaced with the offspring (new population).
        while len(offspring) < len(population):

            # Extract the parents from the population.
            parents = toolbox.select(population, 2)
            parent1, parent2 = toolbox.clone(parents[0]), toolbox.clone(parents[1])

            # Check if there should be mating.
            if random.random() < cxpb:
                
                print("Mated")

                # Mate the 2 parents.
                toolbox.mate(parent1,parent2)

                # Delete the old parents.
                del parent1.fitness.values
                del parent2.fitness.values

            # Check if mutation should be applied to parent1.
            if random.random() < mutpb:
                
                print("Mutated")

                toolbox.mutate(parent1)
                del parent1.fitness.values
            
            # Check if mutation should be applied to parent2.
            if random.random() < mutpb:

                print("Mutated")

                toolbox.mutate(parent2)
                del parent2.fitness.values

            # Add to offspring list
            offspring.append(parent1)
            if len(offspring) < len(population):
                offspring.append(parent2)

        print("Offspring population created")

        # Evaluate crossover/mutated individuals.
        evaluation_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, evaluation_ind)
        for ind, fit in zip(evaluation_ind, fitnesses):
            ind.fitness.values = fit

        if len(hof) > 0:

            worst = max(offspring, key=lambda ind: ind.fitness)
            offspring[offspring.index(worst)] = toolbox.clone(hof[0])

        
        record = stats.compile(offspring)
        logbook.record(gen=gen, nevals=len(evaluation_ind), **record)
        print("==========================================================")
        print(logbook.stream)

        # Save the best indivudal for the current generation.
        best_individuals.append(hof[0].fitness)

        # Replace the population with the offspring.
        population[:] = offspring

    return population, logbook, hof        

if __name__ == "__main__":
    
    pop, log, hof = main()
    print(best_individuals)
    print("\n\n==========================================================\nThe best individual is: %s\nwith fitness: %s\n" \
    "==========================================================" % (hof[0].name, hof[0].fitness))
    vol_tot,vol_vessel,vol_nozzle,vol_pad = hof[0].calc_volume()
    print("The volumes for the different parts of the best individual:\n"
    "==========================================================\n"
    "The volume of vessel: %s m^3\n" \
    "The volume of the nozzle: %s m^3\n" \
    "The volume of the pad: %s m^3\n" \
    "--------------------------------\n" \
    "The total volume: %s m^3\n" % (vol_vessel,vol_nozzle,vol_pad,vol_tot))

    mass_tot = hof[0].calc_mass()
    print("=========================================\n" \
    "The total mass of the best individual is: %s Kg" % mass_tot)
    average_time = np.mean(time_instances_tot)
    print("==========================================================")
    print("The average time of the prediction process was: %s" %average_time)

    gen, avg, min_, max_ = log.select("gen", "avg", "min", "max")

    print("==========================================================")
    print("The difference between the best individual at gen 0 and gen " + str(gen[-1]) + " was: " + str(best_individuals[0].values[0]-best_individuals[-1].values[0]))

    plt.plot(gen, avg, label="average, pop-size=" + str(number_of_individuals),linewidth=2)
    plt.plot(gen, min_, label="minimum, pop-size=" + str(number_of_individuals),linewidth=2)
    plt.plot(gen, max_, label="maximum, pop-size=" + str(number_of_individuals),linewidth=2)
    plt.title("The variation of fitness for a population over %s generations" % gen[-1])
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.legend(loc="center left")

    ax = plt.gca()

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.savefig("Figures/pop_" + str(number_of_individuals) + "_gen_" + str(gen[-1]) + "_mut_" + str(0.01) + "_mate_0.01_1D.svg",format='svg')

    plt.show()