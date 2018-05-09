import random
import copy

class slice:
    def __init__(self, number_of_chrs = 5, len_of_chrs = 1000, mutating_level = 0.3, start_time = 0.0, end_time = 1.0):
        self.number_of_chrs = number_of_chrs
        self.len_of_chrs = len_of_chrs
        self.clones = [[0, 1.0]]
        self.segments = {0 : []}
        for i in range(1, self.number_of_chrs + 1):
            self.segments[0].append([i, 1, self.len_of_chrs, 2, 1, start_time, end_time])
        self.positions = {0 : []}
        self.mutating_level = mutating_level


    # Metod of adding clone inside of slice object.
    # Metod take fraction of new clone (0-1.0), number of parental clone to copy all previous CNA and SNV events
    # and fraction strategy which can be of three type.
    # 'static' - default strategy, is take available fraction from non-abberat cell fraction.
    # 'parental' create new clone as fraction from parent clone which should be specified also.
    # 'weighted' create new clone as independent biological process of adding sampled tissue, all other fractions are recalculated.
    def add_clone(self, new_fraction, parent_clone = 0, fraction_strategy = 'static'):
        if fraction_strategy == 'static':
            self.clones[0][1] -= new_fraction
            self.clones.append([len(self.clones), new_fraction])
            if self.clones[0][1] < 0: print('--- Error: there are no enough fractions ---')
        elif fraction_strategy == 'parental':
            self.clones.append([len(self.clones), self.clones[parent_clone][1] * new_fraction])
            self.clones[parent_clone][1] = self.clones[parent_clone][1] * (1.0 - new_fraction)
            if self.clones[0][1] < 0: print('--- Error: there are no enough fractions ---')
        elif fraction_strategy == 'weighted':
            for i in self.clones:
                self.clones[i][1] = self.clones[i][1] * (1.0 - new_fraction)
            self.clones.append([len(self.clones), new_fraction])
        else:
            print('--- Error: wrong fraction strategy specified ---')

        self.segments[self.clones[-1][0]] = copy.deepcopy(self.segments[self.clones[parent_clone][0]])

        self.positions[self.clones[-1][0]] = copy.deepcopy(self.positions[self.clones[parent_clone][0]])

    # Metod of removing clone inside of slice object.
    # Metod take clone (number) which should be removed (fraction become 0) and fraction strategy which can be of three type.
    # 'static' - default strategy, returns available fraction to non-abberat cell fraction.
    # 'parental' returns fraction to parent clone which should be specified also.
    # 'weighted' returns the fraction's share to all other fractions correspondently it's amound.
    def remove_clone(self, clone, parent_clone = 0, fraction_strategy='static'):
        if fraction_strategy == 'static':
            self.clones[0][1] += self.clones[clone][1]
            self.clones[clone][1] = 0.0
            if self.clones[0][1] > 1: self.clones[0][1] = 1.0
        elif fraction_strategy == 'parental':
            self.clones[parent_clone][1] += self.clones[parent_clone][1] * self.clones[clone][1]
            self.clones[clone][1] = 0.0
            if self.clones[parent_clone][1] > 1: self.clones[parent_clone][1] = 1.0
        elif fraction_strategy == 'weighted':
            for i in self.clones:
                self.clones[i][1] = self.clones[i][1] / (1.0 - self.clones[clone][1])
            self.clones[clone][1] = 0.0
        else:
            print('--- Error: wrong fraction strategy specified ---')

    # Method of duplication in one clone.
    # It take clone number where duplication should happen, chromosome, start and end position of duplicate segment and allele.
    # By default start and end position of duplication corresponds to start and end of chromosome.
    # Allele can be 'A' and 'B'. 'A' is default option.
    # Method creates corresponded segments in all clones and change the segments in the target clone.
    def duplication(self, clone, chromosome, start_position = 1, end_position = -1, allele = 'A', start_time = 0.0, end_time = 1.0):
        if end_position == -1:
            end_position = self.len_of_chrs
        start_segment = 0
        end_segment = 0
        for clone_number in self.segments.keys():
            for num, segment in enumerate(self.segments[clone_number]):
                if segment[0] != chromosome: continue
                if segment[2] == start_position: start_position += 1
                if segment[1] == start_position:
                    start_segment = num
                    break
                if segment[1] < start_position < segment[2]:
                    start_segment = num + 1
                    self.segments[clone_number].insert(num + 1, copy.deepcopy(segment))
                    self.segments[clone_number][num][2] = start_position - 1
                    self.segments[clone_number][num + 1][1] = start_position
                    break
            for num, segment in enumerate(self.segments[clone_number]):
                if segment[0] != chromosome: continue
                if segment[2] == end_position - 1: end_position -= 1
                if segment[2] == end_position:
                    end_segment = num
                    break
                if segment[1] < end_position < segment[2]:
                    end_segment = num
                    self.segments[clone_number].insert(num, copy.deepcopy(segment))
                    self.segments[clone_number][num][2] = end_position
                    self.segments[clone_number][num + 1][1] = end_position + 1
                    break

        for i in range(start_segment, end_segment + 1):
            if self.segments[clone][i][3] == 0: print('--- Error: there is not any chromosome to duplication ---')
            self.segments[clone][i][3] += 1
            self.segments[clone][i][5] = start_time
            self.segments[clone][i][6] = end_time
        if allele == 'B':
            for i in range(start_segment, end_segment+ 1):
                if self.segments[clone][i][3] == 0: print('--- Error: there is not any B chromosome to duplication ---')
                self.segments[clone][i][4] += 1

    # Metod of deletion in one clone.
    # It take clone number where deletion should happen, chromosome, start and end position of delete segment and allele.
    # By default start and end position of deletion corresponds to start and end of chromosome.
    # Allele can be 'A' and 'B'. 'B' is default option.
    # Method creates corresponded segments in all clones and change the segments in the target clone.
    def deletion(self, clone, chromosome, start_position = 1, end_position = -1, allele = 'B', start_time = 0.0, end_time = 1.0):
        if end_position == -1:
            end_position = self.len_of_chrs
        start_segment = 0
        end_segment = 0
        for clone_number in self.segments.keys():
            for num, segment in enumerate(self.segments[clone_number]):
                if segment[0] != chromosome: continue
                if segment[2] == start_position: start_position += 1
                if segment[1] == start_position:
                    start_segment = num
                    break
                if segment[1] < start_position < segment[2]:
                    start_segment = num + 1
                    self.segments[clone_number].insert(num + 1, copy.deepcopy(segment))
                    self.segments[clone_number][num][2] = start_position - 1
                    self.segments[clone_number][num + 1][1] = start_position
                    break
            for num, segment in enumerate(self.segments[clone_number]):
                if segment[0] != chromosome: continue
                if segment[2] == end_position - 1: end_position -= 1
                if segment[2] == end_position:
                    end_segment = num
                    break
                if segment[1] <= end_position <= segment[2]:
                    end_segment = num
                    self.segments[clone_number].insert(num, copy.deepcopy(segment))
                    self.segments[clone_number][num][2] = end_position
                    self.segments[clone_number][num + 1][1] = end_position + 1
                    break

        for i in range(start_segment, end_segment + 1):
            self.segments[clone][i][3] -= 1
            self.segments[clone][i][5] = start_time
            self.segments[clone][i][6] = end_time
            if self.segments[clone][i][3] < 0: print('--- Error: there is not chromosome to deletion ---')
        if allele == 'B':
            for i in range(start_segment, end_segment +1):
                self.segments[clone][i][4] -= 1
                if self.segments[clone][i][4] < 0: print('--- Error: there is not B chromosome to deletion ---')


    # Method of mutating DNA in all clones. Each clone is mutated independently.Mutation happeds accordin time preiod and mutation level.
    # It take value of start and end time point of current period from tissue developed history.
    # Method creates table of number snv according to spesified mutation level.
    def mutating(self, start_time, end_time):
        mutation_probability = (end_time - start_time) * self.mutating_level
        for clone_name in self.positions.keys():
            for chromosome in range(1, self.number_of_chrs + 1):
                for position in range(1, self.len_of_chrs +1):
                    if random.random() < mutation_probability:
                        self.positions[clone_name].append([chromosome, position, start_time, end_time, random.choice('AB')])
            self.positions[clone_name] = sorted(self.positions[clone_name], key = lambda x: (x[0], x[1]))

    # Method to create data of segments and it's ploidy and baf
    def createSegmentsData(self):
        segment_data = []
        for CNA_number in range(len(self.segments[0])):
            mean_ploidy = 0
            mean_baf = 0
            for clone in self.clones:
                mean_ploidy += (self.segments[clone[0]][CNA_number][3] * clone[1])
                mean_baf += (self.segments[clone[0]][CNA_number][4] * clone[1])
            segment_data.append([self.segments[0][CNA_number][0], self.segments[0][CNA_number][1], self.segments[0][CNA_number][2], mean_ploidy - mean_baf, mean_baf])
        return segment_data

    '''
    # Metod to create data of SNV position by timepoint
    def createPositionData(self):
        position_events = []
        for clone in self.clones:
            for position in self.positions[clone[0]]:
                fraction_by_time = random.uniform(1 - position[3], 1 - position[2])
                ploidy = 0.0
                for segment in self.segments[clone[0]]:
                    if position[0] != segment[0]: continue
                    if segment[1] <= position[1] <= segment[2]:
                        if position[4] == 'A': cna_component = segment[3]-segment[4]
                        if position[4] == 'B': cna_component = segment[4]
                        if position[3] <= segment[6]: cna_component = 1.0/segment[3]
                        fraction_by_clone = fraction_by_time * cna_component * clone[1]
                        ploidy = segment[3]
                        break
                position_events.append([position[0], position[1], fraction_by_clone, ploidy, clone[1]])
        position_events = sorted(position_events, key = lambda x: (x[0], x[1]))
        position_data = [position_events[0]]
        position_data[0][2] = 0.0
        position_data[0][3] = 0.0
        for position in position_events:
            if position[0] != position_data[-1][0]:
                position_data.append(position)
                continue
            if position[1] != position_data[-1][1]:
                position_data.append(position)
                continue
            position_data[-1][2] += position[2]
            position_data[-1][3] += position[3] * position[4]
        return position_data
    '''

    # Metod to create data of SNV positionby analytically rules.
    # It prodused list of elements, each of it is of list containing chromosome, position,
    def createPositionData(self):
        position_events = []
        for clone in self.clones:
            for position_number in range(len(self.positions[clone[0]])):
                sum_ploidy = 0.0
                for segment_number in range(len(self.segments[clone[0]])):
                    if self.positions[clone[0]][position_number][0] != self.segments[clone[0]][segment_number][0]: continue
                    if self.segments[clone[0]][segment_number][1] <= self.positions[clone[0]][position_number][1] <= self.segments[clone[0]][segment_number][2]:
                        for cl in self.clones:
                            sum_ploidy += self.segments[cl[0]][segment_number][3] * cl[1]
                        if self.positions[clone[0]][position_number][4] == 'A':
                            ploidy_by_cna = self.segments[clone[0]][segment_number][3] - self.segments[cl[0]][segment_number][4]
                        if self.positions[clone[0]][position_number][4] == 'B':
                            ploidy_by_cna = self.segments[clone[0]][segment_number][4]
                fraction_by_clone = ploidy_by_cna * clone[1]
                position_events.append([self.positions[clone[0]][position_number][0], self.positions[clone[0]][position_number][1], fraction_by_clone, sum_ploidy])
        position_events = sorted(position_events, key=lambda x: (x[0], x[1]))
        print(position_events)
        position_data = [position_events[0]]
        position_data[0][2] = 0.0
        position_data[0][3] = 0.0
        for position in position_events:
            if position[0] != position_data[-1][0]:
                position_data.append(position)
                continue
            if position[1] != position_data[-1][1]:
                position_data.append(position)
                continue
            position_data[-1][2] += position[2]
        position_data.pop(0)
        return position_data

