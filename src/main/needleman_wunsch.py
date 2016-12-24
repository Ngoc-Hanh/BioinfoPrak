from prakt.nw import NeedlemanWunschBase
import random


@NeedlemanWunschBase.register
class NeedlemanWunsch(NeedlemanWunschBase):
    """Needleman Wunsch algorithm"""
    __segA = ""
    __segB = ""
    __lengthA = 0
    __lengthB= 0
    __similarityMatrix ={}
    __tracebackMatrix = {}
    __gapCost = 0
    __complete = False
    __substitutionMatrix = {}
    __alignA = ""
    __alignB = ""
    __alignment = ""

    def initialize(self):

        self.__similarityMatrix[0,0] = 0
        for i in range(1,self.__lengthA+1):
            self.__similarityMatrix[i, 0] = self.__similarityMatrix[i-1,0]  + self.__gapCost
            self.__tracebackMatrix[i,0] = ["up"]

        for j in range(1,self.__lengthB+1):
            self.__similarityMatrix[0,j] = self.__similarityMatrix[0, j-1]  + self.__gapCost
            self.__tracebackMatrix[0, j] = ["left"]

        print("done initialization")
        return

    def fillMatrix(self):
        for i in range(0,self.__lengthA):
            for j in range (0,self.__lengthB):
                currentPosition = (i +1, j+1)
                keyA = self.__segA[i]
                keyB = self.__segB[j]
                subtitutionKey = (keyA, keyB)
                extendAKey = (keyA, '*')
                extendBKey = ('*', keyB)

                # fill similarity matrix

                valueLeft = self.__similarityMatrix[i+1,j]  + int(self.__substitutionMatrix[extendBKey])
                valueUp =   self.__similarityMatrix[i,j+1]  + int(self.__substitutionMatrix[extendAKey])
                valueDiag =   self.__similarityMatrix[i,j]  + int(self.__substitutionMatrix[subtitutionKey])

                value = max (valueLeft, valueUp, valueDiag)

                self.__similarityMatrix [currentPosition] = value

                # fill traceback matrix
                traces = []
                if value == valueLeft:
                    traces.append("left")
                if value == valueUp:
                    traces.append("up")
                if value == valueDiag:
                    traces.append("diag")
                self.__tracebackMatrix [currentPosition] = traces
        print ("done filling matrix")
        return

    def traceback(self,complete_traceback):
        if complete_traceback :
            pass
        else:
            findingAlignment = True
            currentPosition = [self.__lengthA, self.__lengthB]
            traces = []
            while findingAlignment:
                if (currentPosition[0] == 0 and currentPosition[1] == 0):
                    findingAlignment = False
                    break
                else:
                    traces = self.__tracebackMatrix[tuple(currentPosition)]
                    traces = random.sample(traces, len(traces))
                    if traces[0] == "up":
                        self.__alignA = self.__segA[currentPosition[0]-1]
                        self.__alignB = "-"
                        self.__alignment = " "
                        currentPosition[0] -= 1
                    elif traces[0] == "left":
                        self.__alignA = "-"
                        self.__alignB = self.__segB[currentPosition[1]-1]
                        self.__alignment = " "
                        currentPosition[1] -= 1

                    else:
                        self.__alignA = self.__segA[currentPosition[0]-1]
                        self.__alignB = self.__segB[currentPosition[1]-1]
                        if self.__alignA == self.__alignB:
                            self.__alignment = "*"
                        else:
                            self.__alignment = ":"
                        currentPosition[0] -= 1
                        currentPosition[1] -= 1
        print("done traceback")

    def run(self,
            seq1_fasta_fn,
            seq2_fasta_fn,
            subst_matrix_fn,
            cost_gap_open,
            complete_traceback):
            """Document me!"""

            self.__segA = seq1_fasta_fn
            self.__segB = seq2_fasta_fn
            self.__lengthA = len(seq1_fasta_fn)
            self.__lengthB = len(seq2_fasta_fn)
            self.__substitutionMatrix = subst_matrix_fn
            self.__gapCost = cost_gap_open
            self.__complete = complete_traceback

            self.initialize()
            self.fillMatrix()
            self.traceback(complete_traceback)
            print("done")
            score= max(self.__similarityMatrix.values())
            return (score, self.__alignA, self.__alignment, self.__alignB)


# if __name__ == '__main__':
#     # run Needleman-Wunsch with some parameters
#     nw = NeedlemanWunsch()
#     nw.run(
#         "data/sequence1.fa",
#         "data/sequence2.fa",
#         "data/blosum62.txt",
#         5,
#         True)
