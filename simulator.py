import sys
import sliceClass
import interfaceFunctions

arg = sys.argv
number_of_chr, len_of_chr, mutation_level, number_of_cna, mean_DP, noise = interfaceFunctions.readMainDescription(arg[1])
simulation = sliceClass.slice(number_of_chrs = number_of_chr, len_of_chrs = len_of_chr, mutating_level = mutation_level)

segment_data_list = []
position_data_list =[]
cna_step = 0
command_list = interfaceFunctions.readCommandDescription(arg[1], number_of_cna)
for command in command_list:
    start_time = cna_step/number_of_cna
    cna_step += 1
    end_time = cna_step/number_of_cna
    simulation.mutating(start_time= start_time, end_time= end_time)
    if command[0] == 'del':
        simulation.deletion(clone=command[1], chromosome=command[2], start_position=command[3], end_position=command[4], allele=command[5], start_time=start_time, end_time=end_time)
    if command[0] == 'dup':
        simulation.duplication(clone=command[1], chromosome=command[2], start_position=command[3], end_position=command[4], allele=command[5], start_time=start_time, end_time=end_time)
    if command[0] == 'add':
        simulation.add_clone(new_fraction=command[1], parent_clone=command[2], fraction_strategy=command[3])
    if command[0] == 'rem':
        simulation.remove_clone(clone=command[1], parent_clone=command[2], fraction_strategy=command[3])
    if command[0] == 'sam':
        segment_data_list.append([command[1], simulation.createSegmentsData()])
        position_data_list.append([command[1], simulation.createPositionData()])

interfaceFunctions.writeSegments(segment_data_list=segment_data_list)
interfaceFunctions.writeVariants(position_data_list=position_data_list, normal_DP=mean_DP, noise_level=noise)
interfaceFunctions.writeZmatrix(simulation.createClonePositionData(), simulation.clones)


