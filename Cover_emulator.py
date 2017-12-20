#This script run with 3 arguments first is comand '-normal',
# second is lengh of DNA data
# and third is mean resd depth of normal sample.
#Then script create template for data with 3 chromosome  and with noted summery DNA lengh and mean read depth.

import sys

def createNormalTemplate(simulation, DNAlengh, meanDP):
    project_file=open(simulation, 'w')
    project_file.write('# '+simulation+'\n# start\tend\tmeanDP\n')
    project_file.write(str(1)+'\t'+str(DNAlengh)+'\t'+str(meanDP))
    project_file.close()

def test_createNormalTemplate():
    createNormalTemplate('test_createNormalTemplate', 20, 10)
    test_file=open('test_createNormalTemplate', 'r')
    for i in test_file:
        if i[0]!='#':
            result = i.strip().split()
    assert result[0] == '1'
    assert result[1] == '20'
    assert result[2] == '10'

arg=sys.argv

if arg[1]=='-normal':
    createNormalTemplate('current_project', arg[2], arg[3])


