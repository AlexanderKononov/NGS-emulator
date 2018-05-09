import random

# Function take name of file with scenario of run and extract main information
# about samles from heder: number and length of chromosome, mutation level.
# Also it estimates time of tissue developing.
def readMainDescription(file_name):
    number_of_chr = 0
    len_of_chr = 0
    mutation_level = 0.0
    number_of_cna = 0
    mean_DP = 0
    noise = -1
    file_r = open(file_name, 'r')
    for line in file_r:
        if line == '': continue
        #line = line.replace('@', '#').replace('!', '#')
        if line[0] == '#':
            line = line.replace('Num', 'num').replace('Len', 'len').replace('Mut', 'mut').replace('Mea', 'mea').replace('Noi', 'noi')
            if line.find('num') != -1:
                line = line.replace('=', ' ').replace(':', ' ')
                line_p = line.split()
                number_of_chr = int(line_p[-1].strip())
            if line.find('len') != -1:
                line = line.replace('=', ' ')
                line_p = line.split()
                len_of_chr = int(line_p[-1].strip())
            if line.find('mut') != -1:
                line = line.replace('=', ' ')
                line_p = line.split()
                mutation_level = float(line_p[-1].strip().replace(',', '.'))
            if line.find('mea') != -1:
                line = line.replace('=', ' ')
                line_p = line.split()
                mean_DP = int(line_p[-1].strip())
            if line.find('noi') != -1:
                line = line.replace('=', ' ')
                line_p = line.split()
                noise = int(line_p[-1].strip())
            continue
        number_of_cna += 1
    file_r.close()
    if number_of_chr == 0:
        number_of_chr = 5
        print('--- Error: there is not mention about number of chromosomes in scenario! 5 was chosen. ---')
    if len_of_chr == 0:
        len_of_chr = 1000
        print('--- Error: there is not mention about length of chromosomes in scenario! 1000 was chosen. ---')
    if mutation_level == 0.0:
        mutation_level = 0.2
        print('--- Error: there is not mention about mutation level in scenario! 0.2 was chosen. ---')
    if mean_DP == 0:
        mean_DP = 50
        print('--- Error: there is not mention about mean DP in scenario! 50 was chosen. ---')
    if noise == -1:
        noise = 3
        print('--- Error: there is not mention about noise in scenario! 3 was chosen. ---')
    return number_of_chr, len_of_chr, mutation_level, number_of_cna, mean_DP, noise

def readCommandDescription(file_name, number_to_check = 0):
    command_list = []
    discripter = {}
    discripter['del'] = {'clon':1, 'chr':2, 'coo':3, 'stran':4}
    discripter['dup'] = {'clon':1, 'chr':2, 'coo':3, 'stran':4}
    discripter['add'] = {'fr':1, 'par':2, 'strat':3}
    discripter['rem'] = {'clon':1, 'par':2, 'strat':3 }
    discripter['sam'] = {'nam':1}
    file_r = open(file_name, 'r')
    for line in file_r:
        if line == '': continue
        #line = line.replace('@', '#').replace('!', '#')
        if line[0] == '#': continue
        line = line.replace('=', ' ')
        line_p = line.strip().split()
        for command_type in discripter.keys():
            if line_p[0].find(command_type) != -1: line_p[0] = command_type
        command = [line_p[0]]
        for attribute in range(len(discripter[line_p[0]])):
            command.append(0)
        for p in range(1, len(line_p)-1):
            for attribute in discripter[line_p[0]].keys():
                if line_p[p].find(attribute) != -1:
                    command[discripter[line_p[0]][attribute]] = line_p[p + 1]
                    break
        command_list.append(command)
    for command in command_list:
        if command[0] == 'del' or command[0] == 'dup':
            coordinates = command[3].strip().split('-')
            command[3] = coordinates[0]
            command.insert(4, coordinates[1])
            for value in range(1, len(command) - 1):
                command[value] = int(command[value])
            continue
        if command[0] == 'add' or command[0] == 'rem':
            command[1] = float(command[1].replace(',', '.'))
            command[2] = int(command[2])
            if command[0] == 'rem': command[1] = int(command[1])
            continue
    if len(command_list) == 0:
        print('--- Error: There are not any commands in scenario! ---')
    if len(command_list) != number_to_check:
        print('--- Error: Mistakes during reading commands from scenario! ---')
    return command_list

# Method to write segment table with cna data for Canopy analysis
def writeSegments(segment_data_list):
    f_w = open('segments.csv', 'w')
    f_w.write('file\tsample\tchr\tstartpos\tendpos\tnMajor\tnMinor\n')
    for segment_data in segment_data_list:
        for segment in segment_data[1]:
            f_w.write(segment_data[0] + '\t' + segment_data[0])
            for colomn in segment: f_w.write('\t' + str(colomn))
            f_w.write('\n')
    f_w.close()

# Method to write variant table with snv data for Canopy analysis
def writeVariants(position_data_list, normal_DP = 50, noise_level = 3):
    f_w = open('variants.csv', 'w')
    f_w.write('chr\tpos\t.\t.\tID')
    for sample in position_data_list:
        f_w.write('\t' + sample[0] + '.AD\t' + sample[0] + '.DP')
    f_w.write('\n')
    for coordinate in range(len(position_data_list[0][1])):
        f_w.write(str(position_data_list[0][1][coordinate][0]) + '\t' + str(position_data_list[0][1][coordinate][1]) + '\t.\t.\tID1111')
        for sample in position_data_list:
            DP = int(random.normalvariate(sample[1][coordinate][3] * normal_DP, noise_level))
            AF = int(random.normalvariate(sample[1][coordinate][2] * normal_DP, noise_level))
            if AF < 0: AF = 1
            f_w.write('\t' + str(DP - AF) + ',' + str(AF) + '\t' + str(DP))
        f_w.write('\n')
    f_w.close()