import sys
from ArgumentParser import SequenceAlignmentParser
import SubstitutionMatrix
import SequenceLists
import needleman_wunsch
import gotoh
import getopt

def main(argv):

    # Manage options
    parser = SequenceAlignmentParser(sys.argv[1:])
    inputFile = parser.get_input_file()
    outputFile = parser.get_output_file()
    algo = parser.get_algo()
    gapCostOpen = parser.get_gap_cost_open()
    gapCostExtend = parser.get_gap_cost_extend()
    substitutionMatrixType = parser.get_substitution_matrix()
    completeTraceback = parser.is_complete_traceback()

    # Load sequences
    path = sys.argv[0]
    dataDirectory = path.find("src")
    directory = path[0:dataDirectory - 1] + "/data/"
    filename = directory + "testFile.txt"

    sequences = SequenceLists.SequenceLists(filename).values()
    seqs = list(sequences)

   # TODO: remember id
    if (len(seqs) == 0):
        print("not enough sequences")
        sys.exit(2)
    seqA = seqs[0]
    seqB = seqs[1]


    # Load substitution matrix
    substitutionMatrix = SubstitutionMatrix.LoadSubstitutionMatix(substitutionMatrixType)

    # run Algorithms
    print("--------------------------")
    if algo == "gotoh":
        parser.print_info()
        gotohAlign = gotoh.Gotoh()
        print("--------------------------")

        result = gotohAlign.run(seqA, seqB, substitutionMatrix, gapCostOpen, gapCostExtend, False)

    else:
        parser.print_info()
        seqAlign = needleman_wunsch.NeedlemanWunsch()
        print("--------------------------")

        result = seqAlign.run(seqA, seqB, substitutionMatrix, gapCostOpen, False)

    print("score = ", result[0])

if __name__ == '__main__':
    main(sys.argv);