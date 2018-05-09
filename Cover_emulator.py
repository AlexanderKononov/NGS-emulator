# This script contine functions, class and code for run itself (executable code).
# It can get command and additional options, suitable to this command, as arguments after script name.
# Now 5 commands are available: '-normal', '-addclone', '-duplication' '-deletion' and '-seq'

import sys
import random

# Command '-normal' is first command by which inner project files are created.
# Originally this files have description for normal sample without any abberation.
# This script run with command '-normal' by 2-3 arguments:
# first is command '-normal',
# second is length of DNA data per chromosome,
# third is optional argument and notes a number of chromosomes (default is 23).

def createNormalSample(simulation, DNAlengh, number_of_chromosome = 23):
    with open(simulation + '_CNA', 'w') as simulation_CNA_file, open(simulation + '_BAF', 'w') as simulation_BAF_file:
        simulation_CNA_file.write('# ' + simulation + ' 1 1.0\n# chr\tstart\tend\tnormalPloidy\n')
        simulation_BAF_file.write('# ' + simulation + ' 1 1.0\n# chr\tstart\tend\tnormalPloidy\n')
        for i in range(1, int(number_of_chromosome)+1):
            simulation_CNA_file.write(str(i) + '\t' + str(1) + '\t' + str(DNAlengh) + '\t2\n')
        for i in range(1, int(number_of_chromosome)+1):
            simulation_BAF_file.write(str(i) + '\t' + str(1) + '\t' + str(DNAlengh) + '\t1\n')
    with open(simulation + '_phylogeny', 'w') as simulation_phylogeny_file:
        simulation_phylogeny_file.write('# ' + simulation + ' 1 1.0\n# clone_heir\tclone_parent\n')


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

# Command '-seq' is last command by which data from inner project files is transformed to output files.
# This script run with command '-seq' by 3 arguments:
# first is command '-seq',
# second is mean read depth of haploid sample and reflects number of read from one copy of DNA matrix.
# and third is level of noise.
# Then script creates data file by template created before.

def sequence_sample(simulation, meanDP, noise):
    container = DataContainer()
    container.read_data(simulation)
    container.read_phylogeny(simulation)
    with open(simulation + '_CNA.data.txt', 'w') as output_CNA, open(simulation + '_BAF.data.txt', 'w') as output_BAF, open(simulation + '_SNV.data.txt', 'w') as output_SNV:
        for i in range(min(len(container.data_CNA), len(container.data_BAF))):
            segment_CNA = container.data_CNA[i]
            segment_BAF = container.data_BAF[i]
            total_ploidy = 0.0
            total_alleles = 0.0
            segment_SNV = []
            levels_of_SNV = {}
            for j in range(len(container.fractions)):
                total_ploidy += segment_CNA[3 + j] * container.fractions[j]
                total_alleles += segment_BAF[3 + j] * container.fractions[j]
            for key in container.phylogeny.keys()
                levels_of_SNV[key] =
            for j in range(0, len(container.fraction), -1):
                levels_of_SNV[j] = 0.5 * container.fractions[j] + levels_of_SNV[container.phylogeny[j]
                levels_of_SNV[j] = container.fractions[j]

            total_baf=total_alleles/total_ploidy
            for position in range(segment_CNA[1], segment_CNA[2]):
                position_DP = int(random.normalvariate(int(total_ploidy * int(meanDP)), int(noise)))
                output_CNA.write(str(segment_CNA[0]) + '\t' + str(position*1000) + '\t' + str(position_DP) + '\t1\n')
                if position%5 == 0:
                    baf = total_baf + (1-2*total_baf)*random.randint(0,1)
                    position_AD = min(int(random.normalvariate((int(total_ploidy * baf * int(meanDP))), int(noise))), position_DP)
                    output_BAF.write(str(segment_BAF[0]) + '\t' + str(position*1000) + '\t' + str(position_AD) + '\t' + str(position_DP) + '\n')



#This is a description of main data object - data container.
#The containers are used as temporal data format by which the data can be flexibly changed and then can back to project files.
#In the container saved all information about ploidy of each sample's segment from current project file.
#Then the container is changed by metods descriped in the class.
#Then by one of this metods changed data from container is rewrited in current project file.

class DataContainer:
    def __init__(self):
        self.name_of_simulation = '-empty container-'
        self.clonality = 0
        self.fractions = []
        self.data_CNA = []
        self.data_BAF = []
        self.phylogeny ={}

    #Metod of saving data to container from project files.
    def read_data(self, simulation):
        with open(simulation + '_CNA', 'r') as simulation_CNA_file, open(simulation + '_BAF', 'r') as simulation_BAF_file:
            string = simulation_CNA_file.readline().strip().split()
            self.name_of_simulation = string[1]
            self.clonality = int(string[2])
            for i in range(3, len(string)):
                self.fractions.append(float(string[i]))
            for line in simulation_CNA_file:
                if line[0]=='#': continue
                self.data_CNA.append([int(x) for x in line.strip().split()])
            for line in simulation_BAF_file:
                if line[0]=='#': continue
                self.data_BAF.append([int(x) for x in line.strip().split()])


    #Metod of writing data to project files from container.
    def write_data(self, simulation):
        with open(simulation + '_CNA', 'w') as simulation_CNA_file, open(simulation + '_BAF', 'w') as simulation_BAF_file:
            simulation_CNA_file.write('# ' + simulation + ' ' + str(self.clonality))
            simulation_BAF_file.write('# ' + simulation + ' ' + str(self.clonality))
            for i in self.fractions:
                simulation_CNA_file.write(' ' + str(i))
                simulation_BAF_file.write(' ' + str(i))
            simulation_CNA_file.write('\n# chr\tstart\tend\tnormalPloidy')
            simulation_BAF_file.write('\n# chr\tstart\tend\tnormalPloidy')
            for i in range(len(self.fractions)-1):
                simulation_CNA_file.write('\tclone-'+str(i+1))
                simulation_BAF_file.write('\tclone-' + str(i + 1))
            simulation_CNA_file.write('\n')
            simulation_BAF_file.write('\n')
            for i in self.data_CNA:
                for j in range(len(i)):
                    simulation_CNA_file.write(str(i[j]))
                    if j == len(i)-1:
                        simulation_CNA_file.write('\n')
                        continue
                    simulation_CNA_file.write('\t')
            for i in self.data_BAF:
                for j in range(len(i)):
                    simulation_BAF_file.write(str(i[j]))
                    if j == len(i)-1:
                        simulation_BAF_file.write('\n')
                        continue
                    simulation_BAF_file.write('\t')

    #Metod of adding clone inside of container.
    def add_clone(self, new_fraction, parent):
        self.clonality += 1
        self.fractions.append(new_fraction)
        self.fractions[parent] -= new_fraction
        for i in self.data_CNA:
            i.append(i[3])
        for i in self.data_BAF:
            i.append(i[3])
        self.phylogeny[self.clonality] = parent

    # Metod of duplication of one chromosome in one clone inside of container.
    def duplication(self, clone, chromosome, allele = 'majore'):
        for i in self.data_CNA:
            if i[0] == chromosome: i[3+clone] += 1
        if allele == 'minor':
            for i in self.data_BAF:
                if i[0] == chromosome: i[3 + clone] += 1

    # Metod of deletion of one chromosome in one clone inside of container.
    def deletion(self, clone, chromosome, allele = 'minor'):
        for i in self.data_CNA:
            if i[0] == chromosome: i[3+clone] -= 1
        if allele == 'minor':
            for i in self.data_BAF:
                if i[0] == chromosome: i[3 + clone] -= 1

    #Method to transfer pylogenetic data of clone to dictionary in the class
    def read_phylogeny(self, simulation):
        with open(simulation + '_phylogeny', 'r') as simulation_phylogeny_file:
            string = simulation_phylogeny_file.readline().strip().split()
            if self.name_of_simulation != string[1]:
                print('phylogeny_error')
            if self.clonality != int(string[2]):
                print(str(self.clonality))
                print(str(string[2]))
                print('clnality_error')
            for line in simulation_phylogeny_file:
                if line[0]=='#': continue
                p_line = line.strip().split()
                self.phylogeny[int(p_line[0])] = int(p_line[1])

    #Method of writing phylogenetic dictionary to project file
    def write_phylogeny(self, simulation):
        with open (simulation + '_phylogeny', 'w') as simulation_phylogeny_file:
            simulation_phylogeny_file.write('# ' + simulation + ' 1 1.0\n# clone_heir\tclone_parent\n')
            for key in self.phylogeny.keys():
                simulation_phylogeny_file.write(str(key) + '\t' + str(self.phylogeny[key]) + '\n')

# This is function of adding clone.
# It is used with 2 arguments: name of simulation and new clone fraction (float value).
# Sum of all clone fraction must be less then 1.0 .

def add_clone(simulation, new_fraction, parent = 0):
    new_container = DataContainer()
    new_container.read_data(simulation)
    new_container.read_phylogeny(simulation)
    new_container.add_clone(float(new_fraction), int(parent))
    new_container.write_data(simulation)
    new_container.write_phylogeny(simulation)

# This is function of duplication.
# It is used with 3 arguments: name of simulation, clone of changing (int value) and chromosome of changing (int value).

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

# Below code is executable code.
# The part contains API to detect conditions (command and additional options) from arguments during a run.

arg = sys.argv
if arg[1] == '-normal':
    if len(arg) > 3:
        createNormalSample('current_project', arg[2], arg[3])
    else:
        createNormalSample('current_project', arg[2])

if arg[1] == '-seq':
    sequence_sample('current_project', arg[2], arg[3])

if arg[1] == '-addclone':
    if len(arg) > 3:
        add_clone('current_project', arg[2], arg[3])
    else:
        add_clone('current_project', arg[2])

if arg[1] == '-duplication':
    if len(arg) > 4:
        duplication('current_project', arg[2], arg[3], arg[4])
    else:
        duplication('current_project', arg[2], arg[3])

if arg[1] == '-deletion':
    if len(arg) > 4:
        deletion('current_project', arg[2], arg[3], arg[4])
    else:
        deletion('current_project', arg[2], arg[3])

