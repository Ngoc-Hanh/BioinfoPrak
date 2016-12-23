import sys
from ArgumentParser import SequenceAlignmentParser
import SubstitutionMatrix
import SequenceLists
import needleman_wunsch
import gotoh

def main(argv):
    #i = len(argv)
    #print("this is the src function with %d arguments: %s" %(i, argv[0]))

    # Manage options
    parser = SequenceAlignmentParser(argv)
    inputFile = parser.get_input_file()
    outputFile = parser.get_output_file()
    algo = parser.get_algo()
    gapCostOpen = parser.get_gap_cost_open()
    gapCostExtend = parser.get_gap_cost_extend()
    substitutionMatrixType = parser.get_substitution_matrix()
    completeTraceback = parser.is_complete_traceback()

    # Load sequences
    sequences = SequenceLists.SequenceLists("E:/PythonWorkspace/tmp/testFile.txt").values()
    seqs = list(sequences)

   # TODO: remember id

    seqA = seqs[0]
    seqB = seqs[1]


    # Load substitution matrix
    substitutionMatrix = SubstitutionMatrix.LoadSubstitutionMatix(substitutionMatrixType)

    # run Algorithms
    print("Running ", end="")
    if algo.find("gotoh"):
        parser.print_info()
        algo = gotoh.Gotoh
        #algo.run(seqA, seqB, substitutionMatrix, gapCostOpen, gapCostExtend, completeTraceback)
        pass
    else:
        parser.print_info()
        seqAlign = needleman_wunsch.NeedlemanWunsch
        result = seqAlign.run(seqA, seqB, substitutionMatrix, gapCostOpen, completeTraceback)
        print("score = ", result[0])


if __name__ == '__main__':
    main(sys.argv);