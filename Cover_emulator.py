# This script take command by first argument.
# Now 2 commands are available: '-normal' and '-seq'

import sys
import random

# This script run with command '-normal' by 3 arguments:
# first is command '-normal',
# second is length of DNA data

# Then script creates template for data with 3 chromosomes and with noted summery DNA length and mean read depth.

def createNormalSample(simulation, DNAlengh):
    simulation_file = open(simulation, 'w')
    simulation_file.write('# ' + simulation + ' 1 1.0\n# chr\tstart\tend\tnormalPloidy\n')
    for i in range(1, 4):
        simulation_file.write(str(i) + '\t' + str(1) + '\t' + str(DNAlengh) + '\t2\n')
    simulation_file.close()


def test_create_normal_sample():
    createNormalSample('test_createNormalSample', 20)
    test_file = open('test_createNormalSample', 'r')
    for i in test_file:
        if i[0] != '#':
            result = i.strip().split()
    assert result[1] == '1'
    assert result[2] == '20'
    assert result[3] == '2'
    test_file.close()


# This script run with command '-seq' by 3 arguments:
# first is command '-seq',
# second is mean read depth of normal sample.
# and third is level of noise.
# Then script creates data file by template created before and.

def sequence_sample(simulation, meanDP, noise):
    container = DataContainer()
    container.read_data(simulation)
    output=open(simulation + '.data.txt', 'w')
    for segment in container.data:
        for position in range(segment[1], segment[2]):
            total_ploidy = 0.0
            for i in range(len(container.fractions)):
                total_ploidy += segment[3+i]*container.fractions[i]
            position_DP = int(random.normalvariate(int(total_ploidy * int(meanDP)), int(noise)))
            output.write(str(segment[0]) + '\t' + str(position) + '\t' + str(position_DP) + '\n')

#This is a description of main data object - data container.
#In the container saved all information about ploidy of each sample's segment from current project file.
#Then the container is changed by metods descriped in the class.
#Then by one of this metods changed data from container is rewrited in current project file.

class DataContainer:
    def __init__(self):
        self.name_of_simulation='-empty container-'
        self.clonality=0
        self.fractions=[]
        self.data=[]

#Metod of saving data to container from file.
    def read_data(self, simulation):
        simulation_file = open(simulation, 'r')
        string = simulation_file.readline().strip().split()
        self.name_of_simulation = string[1]
        self.clonality = int(string[2])
        for i in range(3, len(string)):
            self.fractions.append(float(string[i]))
        for line in simulation_file:
            if line[0]=='#':
                continue
            self.data.append([int(x) for x in line.strip().split()])
        simulation_file.close()

#Metod of writing data to file from container.
    def write_data(self, simulation):
        simulation_file = open(simulation, 'w')
        simulation_file.write('# ' + simulation + ' ' + str(self.clonality))
        for i in self.fractions:
            simulation_file.write(' '+str(i))
        simulation_file.write('\n# chr\tstart\tend\tnormalPloidy')
        for i in range(len(self.fractions)-1):
            simulation_file.write('\tclone-'+str(i+1))
        simulation_file.write('\n')
        for i in self.data:
            for j in range(len(i)):
                simulation_file.write(str(i[j]))
                if j == len(i)-1:
                    simulation_file.write('\n')
                    continue
                simulation_file.write('\t')
        simulation_file.close()

#Metod of adding clone inside of container.
    def add_clone(self, new_fraction):
        self.clonality += 1
        self.fractions.append(new_fraction)
        self.fractions[0] -= new_fraction
        for i in self.data:
            i.append(i[3])

# Metod of duplication of one chromosome in one clone inside of container.
    def duplication(self, clone, chromosome):
        for i in self.data:
            if i[0] == chromosome: i[3+clone] += 1

# Metod of deletion of one chromosome in one clone inside of container.
    def deletion(self, clone, chromosome):
        for i in self.data:
            if i[0] == chromosome: i[3+clone] -= 1

# This is function of adding clone.
#It is used with 2 arguments: name of simulation and new clone fraction (float value).

def add_clone(simulation, new_fraction):
    new_container = DataContainer()
    new_container.read_data(simulation)
    new_container.add_clone(float(new_fraction))
    new_container.write_data(simulation)

# This is function of duplication.
#It is used with 3 arguments: name of simulation, clone of changing (int value) and chromosome of changing (int value).

def duplication(simulation, clone, chromosome):
    new_container = DataContainer()
    new_container.read_data(simulation)
    new_container.duplication(int(clone), int(chromosome))
    new_container.write_data(simulation)

# This is function of deletion.
# It is used with 3 arguments: name of simulation, clone of changing (int value) and chromosome of changing (int value).

def deletion(simulation, clone, chromosome):
    new_container = DataContainer()
    new_container.read_data(simulation)
    new_container.deletion(int(clone), int(chromosome))
    new_container.write_data(simulation)

arg = sys.argv
print(arg)
if arg[1] == '-normal':
    createNormalSample('current_project', arg[2])

if arg[1] == '-seq':
    sequence_sample('current_project', arg[2], arg[3])

if arg[1] == '-addclone':
    add_clone('current_project', arg[2])

if arg[1] == '-duplication':
    duplication('current_project', arg[2], arg[3])

if arg[1] == '-deletion':
    deletion('current_project', arg[2], arg[3])

