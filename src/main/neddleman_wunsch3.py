from prakt.nw3 import NeedlemanWunsch3Base
from helper.ArgumentParser import SequenceAlignmentParser
from InOutManager.InputManager import InputManager
from InOutManager.OutputManager import FastaOutputManager
from helper.MultiAlignmentResult import MultiAlignmentResult
import sys


@NeedlemanWunsch3Base.register
class NeedlemanWunsch3(NeedlemanWunsch3Base):
    """Needleman Wunsch algorithm"""
    __segA = ""
    __segB = ""
    __segC = ""
    __lengthA = 0
    __lengthB= 0
    __lengthC = 0
    __substitution_matrix = {}
    __similarity_matrix = {}
    __traceback_matrix = {}
    __gap_cost = 0
    __complete = False
    __traceback_list = []
    __current_traces = []

    def sum_of_pair(self, a, b, c ):
        result = (int(self.__substitution_matrix[(a,b)]) +
               int(self.__substitution_matrix[(b,c)]) +
               int(self.__substitution_matrix[(a,c)]))
        #print ("({}, {}, {}) : {}".format(a, b, c, result))
        return result

    def initialize(self):
        self.__similarity_matrix = {}
        self.__traceback_matrix = {}
        self.__similarity_matrix[0, 0, 0] = 0
        for i in range(1, self.__lengthA + 1):
            self.__similarity_matrix[i, 0, 0] = self.__similarity_matrix[i - 1, 0, 0] + 2*self.__gap_cost
            self.__traceback_matrix[i, 0, 0] = ["gap_BC"]
        for j in range(1, self.__lengthB + 1):
            self.__similarity_matrix[0, j, 0] = self.__similarity_matrix[0, j - 1, 0] + 2*self.__gap_cost
            self.__traceback_matrix[0, j,0] = ["gap_AC"]
        for k in range(1, self.__lengthC + 1):
            self.__similarity_matrix[0, 0, k] = self.__similarity_matrix[0, 0, k-1] + 2*self.__gap_cost
            self.__traceback_matrix[0, 0, k] = ["gap_AB"]

        for i in range(1, self.__lengthA + 1):
            for j in range(1, self.__lengthB + 1):
                self.__similarity_matrix[i, j, 0] = self.__similarity_matrix[i-1, j - 1, 0] + self.__gap_cost
                self.__traceback_matrix[i, j, 0] = ["gap_C"]
            for k in range(1, self.__lengthC + 1):
                self.__similarity_matrix[i, 0, k] = self.__similarity_matrix[i-1, 0, k - 1] + self.__gap_cost
                self.__traceback_matrix[i, 0, k] = ["gap_B"]

        for j in range(1, self.__lengthB + 1):
            for k in range(1, self.__lengthC + 1):
                self.__similarity_matrix[0, j, k] = self.__similarity_matrix[0, j-1, k - 1] + self.__gap_cost
                self.__traceback_matrix[0, j, k] = ["gap_AB"]
        print("done initialization")
        return

    def fill_matrices(self):
        for i in range(1, self.__lengthA+1):
            for j in range(1, self.__lengthB+1):
                for k in range(1, self.__lengthC+1):
                    current_position = (i, j, k)
                    keyA = self.__segA[i-1]
                    keyB = self.__segB[j-1]
                    keyC = self.__segC[k-1]

                    # fill similarity matrix
                    gap_A = self.__similarity_matrix[i-1, j, k] + self.sum_of_pair("*", keyB, keyC)
                    gap_B = self.__similarity_matrix[i, j-1, k] + self.sum_of_pair(keyA, "*", keyC)
                    gap_C = self.__similarity_matrix[i, j, k-1] + self.sum_of_pair(keyA, keyB, "*")
                    gap_AB = self.__similarity_matrix[i-1, j, k] + self.sum_of_pair("*", "*", keyC)
                    gap_AC = self.__similarity_matrix[i-1, j, k] + self.sum_of_pair("*", keyB, "*")
                    gap_BC = self.__similarity_matrix[i-1, j, k] + self.sum_of_pair(keyA, "*", "*")
                    no_gap = self.__similarity_matrix[i-1, j, k] + self.sum_of_pair(keyA, keyB, keyC)

                    current_value = max(gap_A, gap_B, gap_C,
                                    gap_AB, gap_BC, gap_AC,
                                    no_gap)
                    self.__similarity_matrix[current_position] = current_value
                    print("{} ".format(self.__similarity_matrix[current_position]), end="")

                    # fill traceback matrix
                    traces = []
                    if current_value == gap_A:
                        traces.append("gap_A")
                    if current_value == gap_B:
                        traces.append("gap_B")
                    if current_value == gap_C:
                        traces.append("gap_C")
                    if current_value == gap_AB:
                        traces.append("gap_AB")
                    if current_value == gap_AC:
                        traces.append("gap_AC")
                    if current_value == gap_BC:
                        traces.append("gap_BC")
                    if current_value == no_gap:
                        traces.append("no_gap")

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
            if self.__lengthA + self.__lengthB +self.__lengthC<= 10000:
                sys.setrecursionlimit(100000)
                self.recursive_traceback(self.__lengthA, self.__lengthB, self.__lengthC, traceback, complete_traceback)
            else:
                print("the sequences are too large for recursive traceback!")
        alignment_list = []
        for t in self.traceback_list:
            alignment = self.expand_traceback(t)
            alignment_list.append(alignment)
        print("done traceback")
        return alignment_list

    def recursive_traceback(self, currentPosA, currentPosB, currentPosC, traceback, complete_traceback):
        done = False
        if done: return
        if currentPosA == 0 and currentPosB == 0:
            self.traceback_list.append(traceback[:])
            if complete_traceback:
                done = True
                return
            if len(self.__current_traces) >0:
                del traceback[:]
                #traceback.clear()
                traceback.append(self.__current_traces.pop())
            print("list size: {}".format(len(self.traceback_list)))

        # elif currentPosA == 0:
        #      traceback.append("gap_A")
        #      currentPosB -=1
        #      #print("{}, {}". format(currentPosA, currentPosB))
        #      print(traceback)
        #      self.recursive_traceback(currentPosA, currentPosB, traceback)
        #
        # elif currentPosB == 0:
        #      traceback.append("up")
        #      currentPosA -= 1
        #      print(traceback)
        #      self.recursive_traceback(currentPosA, currentPosB, traceback)

        else:
            traces = self.__traceback_matrix[(currentPosA, currentPosB, currentPosC)]
            if len(traces) > 1:
                print("pos: " + str((currentPosA,currentPosB, currentPosC)) + ", more:" + str(traces))
                newtrace= traceback[:]
                self.__current_traces.append(newtrace)

            for trace in traces:
                traceback.append(trace)
                if trace == "gap_A":
                    currentPosB -= 1
                    currentPosC -= 1
                elif trace == "gap_B":
                    currentPosA -= 1
                    currentPosC -= 1
                elif trace == "gap_C":
                    currentPosA -= 1
                    currentPosB -= 1
                elif trace == "gap_AB":
                    currentPosC -= 1
                elif trace == "gap_AC":
                    currentPosB -= 1
                elif trace == "gap_BC":
                    currentPosA -= 1
                else:
                    currentPosA -= 1
                    currentPosB -= 1
                    currentPosC -= 1
                print(traceback)
                self.recursive_traceback(currentPosA, currentPosB, currentPosC, traceback, complete_traceback)

    def expand_traceback(self, traceback):
        alignA = ""
        alignB = ""
        alignC = ""
        a = 0
        b = 0
        c = 0
        for trace in traceback:
            if trace == "gap_A":
                alignA = alignA + "-"
                alignB = alignB + self.__segB[b]
                alignC = alignC  +self.__segC[c]
                b += 1
                c += 1
            elif trace == "gap_B":
                alignA = alignA + self.__segA[a]
                alignB = alignB + "-"
                alignC = alignC  +self.__segC[c]
                a += 1
                c += 1

            elif trace == "gap_C":
                alignA = alignA + self.__segA[a]
                alignB = alignB + self.__segB[b]
                alignC = alignC  +"-"
                a += 1
                b += 1
            elif trace == "gap_AB":
                alignA = alignA + "-"
                alignB = alignB + "-"
                alignC = alignC  +self.__segC[c]
                c += 1
            elif trace == "gap_AC":
                alignA = alignA + "-"
                alignB = alignB + self.__segB[b]
                alignC = alignC + "-"
                b += 1
            elif trace == "gap_BC":
                alignA = alignA + self.__segA[a]
                alignB = alignB + "-"
                alignC = alignC  + "-"
                a += 1
            else:
                alignA = alignA + self.__segA[a]
                alignB = alignB + self.__segB[b]
                alignC = alignC  +self.__segC[c]
                a += 1
                b += 1
                c += 1
        return (alignA, alignB, alignC)

    def run_algorithm(self, key1, seq1, key2, seq2, key3, seq3, matrix, gap_cost, complete_traceback):
        self.__segA = seq1
        self.__segB = seq2
        self.__segC =seq3
        self.__lengthA = len(seq1)
        self.__lengthB = len(seq2)
        self.__lengthC = len (seq3)
        self.__gap_cost = gap_cost
        self.__complete = complete_traceback

        # run algorithm
        self.initialize()
        self.fill_matrices()
        score = self.__similarity_matrix[(self.__lengthA, self.__lengthB, self.__lengthC)]

        alignment_list = self.find_tracebacks(complete_traceback)
        print("done")
        result = MultiAlignmentResult()
        result.id = [key1, key2, key3]
        result.score = score
        result.aligned_sequence_list = alignment_list
        return result
    def run(self,
            seq_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            complete_traceback):
        """
        Calculate optimal alignment(s) with Needleman-Wunsch algorithm for three sequences.

        Args:
            seq_fasta_fn: path to fasta file containing sequences
            subst_matrix_fn: path to substitution matrix
            cost_gap_open: cost to open a gap
            complete_traceback: If True, return all optimal alignments. Otherwise choose a random alignment.

        Returns:
            tuple of
            (id_seq1: fasta id of first sequence,
             seq1: first sequence,
             id_seq2: fasta id of second sequence,
             seq2: second sequence,
             id_seq3: fasta id of second sequence,
             seq3: second sequence,
             score: score of optimal alignment,
             [(aln_string_seq1, aln_string_seq2, aln_string_seq3), ...]: list of tuples containing optimal alignments)
        """
        # load inputs
        sequences = {}
        input_manager = InputManager()
        if seq_fasta_fn is None:
            print("No input file! Program will now exit")
            exit(1)
        else:
            sequences = input_manager.load_fasta_sequences(seq_fasta_fn)

        if len(sequences) < 3:
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
        result = self.run_algorithm(keys[0], sequences.get(keys[0]),
                                    keys[1], sequences.get(keys[1]),
                                    keys[2], sequences.get(keys[2]),
                                    self.__substitution_matrix, cost_gap_open, complete_traceback)
        return (keys[0], sequences.get(keys[0]),
                keys[1], sequences.get(keys[1]),
                keys[2], sequences.get(keys[2]),
                result.score,
                result.aligned_sequence_list)

if __name__ == '__main__':
    # run Needleman-Wunsch with some parameters
    nw = NeedlemanWunsch3()
    parser = SequenceAlignmentParser(sys.argv[1:])
    input_options = parser.options
    if input_options is None:
        print('Usage:'
              'needleman_wunsch3.py -i <inputfile1> -n <inputfile2> -m <substitutionMatrix> '
              '                     [-o <outputfile> -g <gapCostOpen> -c]'
              'Optional:'
              '-o    outputfile'
              '-g    linear gap cost'
              '-c    complete traceback')

        exit(1)
    else:
        results = nw.run(input_options.path_to_fasta_sequence_file,
                         input_options.path_to_substitution_matrix,
                         input_options.open_gap_cost,
                         input_options.traceback_complete)
        result = results[7]
        print(results)
        output = FastaOutputManager()
        for i in range(len(results[7])):
            print("-"*10 + "solution {}".format(i) + "-"*10)
            output.display(results[6],
                           results[0], results[7][i][0],
                           results[2], results[7][i][1],
                           results[4],results[7][i][2] )

        if input_options.path_to_output is not None:
            output.save(results[6],
                           results[0], results[7][i][0],
                           results[2], results[7][i][1],
                           results[4],results[7][i][2],
                        input_options.path_to_output)