# This script take command by first argument.
# Now 2 commands are available: '-normal' and '-seq'

import sys
import random

# This script run with command '-normal' by 3 arguments:
# first is command '-normal',
# second is length of DNA data
# and third is mean resd depth of normal sample.
# Then script creates template for data with 3 chromosomes and with noted summery DNA length and mean read depth.

def createNormalSample(simulation, DNAlengh, meanDP):
    project_file = open(simulation, 'w')
    project_file.write('# ' + simulation + '\n# chr\tstart\tend\tmeanDP\n')
    for i in range(1, 4):
        project_file.write(str(i) + '\t' + str(1) + '\t' + str(DNAlengh) + '\t' + str(meanDP) + '\n')
    project_file.close()


def test_createNormalSample():
    createNormalSample('test_createNormalSample', 20, 10)
    test_file = open('test_createNormalSample', 'r')
    for i in test_file:
        if i[0] != '#':
            result = i.strip().split()
    assert result[1] == '1'
    assert result[2] == '20'
    assert result[3] == '10'

# This script run with command '-seq' by 2 arguments:
# first is command '-seq',
# second is level of noise.
# Then script creates data file by template created before and.

def sequenceSample(simulation, noise):
    with open(simulation, 'r') as input, open(simulation + '.data.txt', 'w') as output:
        for line in input:
            if line[0] == '#':
                continue
            position_data = line.strip().split()
            for position in range(int(position_data[1]) + 1, int(position_data[2]) + 1):
                position_DP = int(random.normalvariate(int(position_data[3]), int(noise)))
                output.write(position_data[0] + '\t' + str(position) + '\t' + str(position_DP) + '\n')


arg = sys.argv

if arg[1] == '-normal':
    createNormalSample('current_project', arg[2], arg[3])

if arg[1] == '-seq':
    sequenceSample('current_project', arg[2])
