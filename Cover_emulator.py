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
    simulation_file.write('# ' + simulation + '\n# chr\tstart\tend\tmeanDP\n')
    for i in range(1, 4):
        simulation_file.write(str(i) + '\t' + str(1) + '\t' + str(DNAlengh) + '\t2\n')
    simulation_file.close()


def test_createNormalSample():
    createNormalSample('test_createNormalSample', 20)
    test_file = open('test_createNormalSample', 'r')
    for i in test_file:
        if i[0] != '#':
            result = i.strip().split()
    assert result[1] == '1'
    assert result[2] == '20'
    assert result[3] == '2'

# This script run with command '-seq' by 3 arguments:
# first is command '-seq',
# second is mean read depth of normal sample.
# and third is level of noise.
# Then script creates data file by template created before and.

def sequenceSample(simulation, meanDP, noise):
    with open(simulation, 'r') as input, open(simulation + '.data.txt', 'w') as output:
        for line in input:
            if line[0] == '#':
                continue
            position_data = line.strip().split()
            for position in range(int(position_data[1]) + 1, int(position_data[2]) + 1):
                position_DP = int(random.normalvariate(int(position_data[3])*int(meanDP), int(noise)))
                output.write(position_data[0] + '\t' + str(position) + '\t' + str(position_DP) + '\n')


arg = sys.argv

if arg[1] == '-normal':
    createNormalSample('current_project', arg[2], arg[3])

if arg[1] == '-seq':
    sequenceSample('current_project', arg[2])
