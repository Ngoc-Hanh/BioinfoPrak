from prakt.nw import NeedlemanWunschBase
from helper.ArgumentParser import SequenceAlignmentParser
from InOutManager.InputManager import InputManager
from InOutManager.OutputManager import PairwiseOutputManager
from helper.PairwiseAlignmentResult import PairwiseAlignmentResult
import random
import sys


@NeedlemanWunschBase.register
class NeedlemanWunsch(NeedlemanWunschBase):
    """Needleman Wunsch algorithm"""
    __segA = ""
    __segB = ""
    __lengthA = 0
    __lengthB = 0
    __substitution_matrix = {}
    __similarity_matrix = {}
    __traceback_matrix = {}
    __gap_cost = 0
    __complete = False
    __traceback_list = []
    __current_traces = []
    def initialize(self):
        self.__similarity_matrix = {}
        self.__traceback_matrix = {}
        self.__similarity_matrix[0, 0] = 0
        for i in range(1, self.__lengthA + 1):
            self.__similarity_matrix[i, 0] = self.__similarity_matrix[i - 1, 0] + self.__gap_cost

            self.__traceback_matrix[i, 0] = ["up"]

        for j in range(1, self.__lengthB + 1):
            self.__similarity_matrix[0, j] = self.__similarity_matrix[0, j - 1] + self.__gap_cost
            self.__traceback_matrix[0, j] = ["left"]

        print("done initialization")
        return

    def fill_matrices(self):
        for i in range(1, self.__lengthA+1):
            for j in range(1, self.__lengthB+1):
                current_position = (i, j)
                keyA = self.__segA[i-1]
                keyB = self.__segB[j-1]
                substitution_key = (keyA, keyB)
                extend_A_Key = (keyA, '*')
                extend_B_Key = ('*', keyB)

                # fill similarity matrix
                valueLeft = self.__similarity_matrix[i, j-1] + int(self.__substitution_matrix[extend_B_Key])
                valueUp = self.__similarity_matrix[i-1, j ] + int(self.__substitution_matrix[extend_A_Key])
                valueDiag = self.__similarity_matrix[i-1, j-1] + int(self.__substitution_matrix[substitution_key])

                current_value = max(valueLeft, valueUp, valueDiag)
                self.__similarity_matrix[current_position] = current_value
                print("{} ".format(self.__similarity_matrix[current_position]), end="")

                # fill traceback matrix
                traces = []
                if current_value == valueLeft:
                    traces.append("left")
                if current_value == valueUp:
                    traces.append("up")
                if current_value == valueDiag:
                    traces.append("diag")
                self.__traceback_matrix[current_position] = traces
        print("done filling matrix")
        return

    def find_tracebacks(self, complete_traceback):
        '''
        Calculate tracebacks
        :param complete_traceback: if true, return a list of all possible tracebacks, else a list of a random traceback
        :return: a list of tracebacks
        '''
        traceback = []
        self.traceback_list = []
        self.__current_traces =[]
        if complete_traceback:
            if self.__lengthA + self.__lengthB <= 10000:
                sys.setrecursionlimit(100000)
                self.recursive_traceback(self.__lengthA, self.__lengthB, traceback)
            else:
                print("the sequences are too large for recursive traceback!")
        else:
            findingAlignment = True
            currentPosition = [self.__lengthA, self.__lengthB]
            traces = []
            while findingAlignment:
                if (currentPosition[0] == 0 and currentPosition[1] == 0):
                    findingAlignment = False
                    break
                else:
                    traces = self.__traceback_matrix[tuple(currentPosition)]
                    if len(traces) > 1:
                        print("pos: " + str(currentPosition) + ", more:" + str(traces))
                    traces = random.sample(traces, len(traces))
                    traceback.append(traces[0])
                    if traces[0] == "up":
                        currentPosition[0] -= 1
                    elif traces[0] == "left":
                        currentPosition[1] -= 1
                    else:
                        currentPosition[0] -= 1
                        currentPosition[1] -= 1
                    print(traceback)
            self.traceback_list.append(traceback)
        alignment_list = []
        for t in self.traceback_list:
            alignment = self.expand_traceback(t)
            alignment_list.append(alignment)
        print("done traceback")
        return alignment_list

    def recursive_traceback(self, currentPosA, currentPosB, traceback):
        if currentPosA == 0 and currentPosB == 0:
            self.traceback_list.append(traceback[:])
            if len(self.__current_traces) >0:
                del traceback[:]
                #traceback.clear()
                traceback.append(self.__current_traces.pop())
            print("list size: {}".format(len(self.traceback_list)))

        elif currentPosA == 0:
             traceback.append("left")
             currentPosB -=1
             #print("{}, {}". format(currentPosA, currentPosB))
             print(traceback)
             self.recursive_traceback(currentPosA, currentPosB, traceback)

        elif currentPosB == 0:
             traceback.append("up")
             currentPosA -= 1
             print(traceback)
             self.recursive_traceback(currentPosA, currentPosB, traceback)

        else:
            traces = self.__traceback_matrix[(currentPosA, currentPosB)]
            if len(traces) > 1:
                print("pos: " + str((currentPosA,currentPosB)) + ", more:" + str(traces))
                newtrace= traceback[:]
                self.__current_traces.append(newtrace)

            for trace in traces:
                traceback.append(trace)
                if trace == "up":
                    currentPosA -= 1
                elif trace == "left":
                    currentPosB -= 1
                else:
                    currentPosA -= 1
                    currentPosB -= 1
                print(traceback)
                self.recursive_traceback(currentPosA, currentPosB, traceback)

    def expand_traceback(self, traceback):
        alignA = ""
        alignB = ""
        a = 0
        b = 0
        for trace in traceback:
            if trace == "up":
                alignA = alignA + self.__segA[a]
                alignB = alignB + "-"
                a +=1
            elif trace == "left":
                alignA = alignA + "-"
                alignB = alignB + self.__segB[b]
                b += 1
            else:
                alignA = alignA + self.__segA[a]
                alignB = alignB + self.__segB[b]
                a += 1
                b += 1
        return (alignA, alignB)

    def run_algorithm(self, key1, seq1, key2, seq2, matrix, gap_cost, complete_traceback):
        self.__segA = seq1
        self.__segB = seq2
        self.__lengthA = len(seq1)
        self.__lengthB = len(seq2)
        self.__gap_cost = gap_cost
        self.__complete = complete_traceback

        # run algorithm
        self.initialize()
        self.fill_matrices()
        score = self.__similarity_matrix[(self.__lengthA, self.__lengthB)]

        alignment_list = self.find_tracebacks(complete_traceback)
        print("done")
        result = PairwiseAlignmentResult()
        result.alignment_result_list(key1, key2, score, alignment_list)
        return result

    def run(self,
            seq1_fasta_fn,
            seq2_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm.

        Args:
            seq1_fasta_fn: path to fasta file containing first sequence
            seq2_fasta_fn: path to fasta file containing second sequence
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.

        Returns:
            tuple of
            (id_seq1: fasta id of first sequence,
             seq1: first sequence,
             id_seq2: fasta id of second sequence,
             seq2: second sequence,
             score: score of optimal alignment,
             [(aln_string_seq1, aln_string_seq2), ...]: list of tuples containing optimal alignments)
        """

        # load inputs
        sequences = {}
        input_manager = InputManager()
        if seq1_fasta_fn is None and seq2_fasta_fn is None:
            print("No input file! Program will now exit")
            exit(1)
        else:
            if seq1_fasta_fn == seq2_fasta_fn or seq2_fasta_fn == None:
                print('there is only one input file')
            else:
                sequences = input_manager.load_fasta_sequences(seq2_fasta_fn)
            sequences.update(input_manager.load_fasta_sequences(seq1_fasta_fn))

        if len(sequences) <2:
            print("Not enough sequences. Program will now exit!")
            exit(1)

        # load substitution matrix
        self.__substitution_matrix = input_manager.load_substitution_matrix(subst_matrix_fn)

        if input_manager.not_valid_substitution_matrix:
            print("substitution matrix is not correct!")
            exit(1)
        if input_manager.not_valid_sequences:
            print("at least one of the sequences is not in fasta format")
            exit(1)

        keys = list(sequences.keys())
        result= self.run_algorithm(keys[0], sequences.get(keys[0]),
                                  keys[1], sequences.get(keys[1]),
                                  self.__substitution_matrix, cost_gap_open, complete_traceback)
        return result.to_tuple(sequences.get(keys[0]), sequences.get(keys[1]))

if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch()
    parser = SequenceAlignmentParser(sys.argv[1:])
    input_options = parser.options
    if input_options is None:
        print('Usage:'
              'needleman_wunsch.py -i <inputfile1> -n <inputfile2> -m <substitutionMatrix> '
              '                     [-o <outputfile> -g <gapCostOpen> -c]'
              'Optional:'
              '-o    outputfile'
              '-g    linear gap cost'
              '-c    complete traceback')

        exit(1)
    else:
        results = nw.run(input_options.path_to_fasta_sequence_file,
                         input_options.path_to_sequence_file2,
                         input_options.path_to_substitution_matrix,
                         input_options.open_gap_cost,
                         input_options.traceback_complete)
        result = results[5]
        print(results)
        output = PairwiseOutputManager()
        for i in range(len(results[5])):
            print("-"*10 + "solution {}".format(i) + "-"*10)
            output.display(results[4], results[5][i][0], results[5][i][1])

        if input_options.path_to_output is not None:
            output.save(result, input_options.path_to_output)
