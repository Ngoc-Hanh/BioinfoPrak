import sys
import getopt
from helper.CommandLineOptions import CommandLineOptions


class SequenceAlignmentParser:
    """
    This is a argument parser for all Sequence Alignment Algorithms
    """
    __inputFile = None
    __inputFile2 = None
    __outputFile = None
    __substitutionMatrix = None
    __gapCostOpen = 1
    __gapCostExtend = 1
    __completeTraceback = False
    options = None
    def __init__ (self, arg_list):
        self.options = CommandLineOptions()
        try:
            opts, args = getopt.getopt(arg_list, "hi:n:o:m:g:e:ca:",[]) #, ["ifile=","ifile2="  "ofile=", "substitutionMatrix", "gapCost=", "completeTraceback", "algo="])
            print(len(opts), len(args))
        except getopt.GetoptError:
            print('Arguments: '
                    '-i               inputfile'
                    '-m               substitution matrix'
                    '-n               second input file (not required, if the first one contains multiple sequences'
                    '-o               outputfile'
                    '-g               open gap cost'
                    '-e               extend gap cost')
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-h':
                print('Arguments: '
                      '-i               inputfile'
                      '-m               substitution matrix'
                      '-n               second input file (not required, if the first one contains multiple sequences'
                      '-o               outputfile'
                      '-g               open gap cost'
                      '-e               extend gap cost')
                sys.exit(2)
            elif opt in ("-i", "--ifile"):
                self.__inputFile = arg
            elif opt in ("-n", "--ifile2"):
                self.__inputFile2 = arg
            elif opt in ("-o", "--ofile"):
                self.__outputFile = arg
            elif opt in ("-m", "--substitutionMatrix"):
                    self.__substitutionMatrix = arg
            elif opt in ("-g", "--gapCostOpen"):
                self.__gapCostOpen = int(arg)
            elif opt in ("-e", "--gapCostExtend"):
                self.__gapCostExtend = int(arg)
            elif opt in ("-c", "--completeTraceback"):
                self.__completeTraceback = True

        self.options.setAlignmentOptions(self.__inputFile,
                                        self.__inputFile2,
                                        self.__substitutionMatrix,
                                        self.__gapCostOpen,
                                        self.__gapCostExtend,
                                        self.__completeTraceback,
                                        self.__outputFile)



    def print_info(self):
        print("substitution matrix = ", self.__substitutionMatrix)
        print("open gap cost = ", self.__gapCostOpen)
        print("extend gap cost = ", self.__gapCostExtend)
        print("complete traceback = ", self.__completeTraceback)
