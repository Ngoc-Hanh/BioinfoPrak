import sys
import getopt


class SequenceAlignmentParser:
    """
    This is a argument parser for all Sequence Alignment Algorithms
    """
    __inputFile = ""
    __outputFile = ""
    __substitutionMatrix = "blosum"
    __gapCostOpen = 1
    __gapCostExtend = 1
    __completeTraceback = False
    __algo = "nw"

    def __init__ (self, arg_list):
        try:
            opts, args = getopt.getopt(arg_list, "hi:o:m:g:e:ca:",[]) #, ["ifile=", "ofile=", "substitutionMatrix", "gapCost=", "completeTraceback", "algo="])
            print(len(opts), len(args))
        except getopt.GetoptError:
            print('SequenceAlignment.py -i <inputfile> -o <outputfile> -g <gapCostOpen> -e <gapCostExtend> -m -c -a <algo>')
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-h':
                print('SequenceAlignment.py -i <inputfile> -o <outputfile> -g <gapCost> -e <gapCostExtend> -m -c -a <algo>')
                sys.exit()
            elif opt in ("-i", "--ifile"):
                self.__inputFile = arg
            elif opt in ("-o", "--ofile"):
                self.__outputFile = arg
            elif opt in ("-m", "--substitutionMatrix"):
                if "pam" in arg:
                    self.__substitutionMatrix = "pam"
                else:
                    self.__substitutionMatrix = "blosum"
            elif opt in ("-g", "--gapCostOpen"):
                self.__gapCostOpen = int(arg)
            elif opt in ("-e", "--gapCostExtend"):
                self.__gapCostExtend = int(arg)
            elif opt in ("-c", "--completeTraceback"):
                self.__completeTraceback = True
            elif opt in ("-a", "--algo"):
                if arg == "gotoh":
                    self.__algo = "gotoh"
                else:
                    self.__algo = "nw"

    def get_output_file(self):
        # return path to the output file
        return self.__outputFile

    def get_input_file(self):
        # return path to the output file
        return self.__inputFile

    def get_substitution_matrix(self):
        # return true if pamMatrix option is given, else use Blosum substitution matrix
        return self.__substitutionMatrix

    def get_gap_cost_open(self):
        # return the gap cost, default =1
        return self.__gapCostOpen

    def get_gap_cost_extend(self):
        # return the gap cost, default =1
        return self.__gapCostExtend

    def is_complete_traceback(self):
        # return the option to calculate all the alignments
        return self.__completeTraceback

    def get_algo(self):
        # return algorithm name
        return self.__algo

    def print_info(self):
        print("algo = " , self.__algo)
        print("substitution matrix = ", self.__substitutionMatrix)
        print("open gap cost = ", self.__gapCostOpen)
        print("extend gap cost = ", self.__gapCostExtend)
        print("complete traceback = ", self.__completeTraceback)
