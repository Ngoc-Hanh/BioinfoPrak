from prakt.gt import GotohBase
import random

@GotohBase.register
class Gotoh(GotohBase):
    """Document me!"""
    __segA = ""
    __segB = ""
    __lengthA = 0
    __lengthB = 0
    __similarityMatrixP = {}
    __similarityMatrixD = {}
    __similarityMatrixQ = {}
    __tracebackMatrixP = {}
    __tracebackMatrixD = {}
    __tracebackMatrixQ = {}
    __gapCostOpen = 0
    __gapCostExtend=0
    __complete = False
    __substitutionMatrix = {}
    __alignA = ""
    __alignB = ""
    __alignment = ""

    def initialize(self):
        maxValue = 1000
        # Matrix D
        self.__similarityMatrixD[0, 0] = 0
        for i in range(1, self.__lengthA + 1):
            self.__similarityMatrixD[i, 0] = self.__similarityMatrixD[i - 1, 0] + self.__gapCostOpen
            #self.__tracebackMatrix[i, 0] = ["up"]

        for j in range(1, self.__lengthB + 1):
            self.__similarityMatrixD[0, j] = self.__similarityMatrixD[0, j - 1] + self.__gapCostOpen
            #self.__tracebackMatrix[0, j] = ["right"]

        # Matrix P
        self.__similarityMatrixP[0, 0] = 0

        for j in range(1, self.__lengthB + 1):
            self.__similarityMatrixP[0, j] = maxValue
            #self.__tracebackMatrix[0, j] = ["right"]

        # Matrix Q
        self.__similarityMatrixQ[0, 0] = 0
        for i in range(1, self.__lengthA + 1):
            self.__similarityMatrixQ[i, 0] = maxValue
            #self.__tracebackMatrix[i, 0] = ["up"]

        return

    def fillMatrix(self):
        for i in range(0, self.__lengthA):
            for j in range(0, self.__lengthB):
                diagonalPosition = (i,j)
                upPosition = (i, j+1)
                rightPosition = (i+1, j)
                currentPosition = (i + 1, j + 1)
                keyA = self.__segA[i]
                keyB = self.__segB[j]
                subtitutionKey = (keyA, keyB)

                # fill P
                case1P = self.__similarityMatrixD[upPosition] + self.__gapCostOpen
                case2P =  self.__similarityMatrixP[upPosition] +self.__gapCostExtend
                valueP = max (case1P, case2P)

                # fill Q
                case1Q = self.__similarityMatrixD[rightPosition] + self.__gapCostOpen
                case2Q = self.__similarityMatrixQ[rightPosition] + self.__gapCostExtend
                valueQ = max(case1Q, case2Q)

                # fill D
                valueD = self.__similarityMatrixD[subtitutionKey] + self.__substitutionMatrix(subtitutionKey)
                value = max(valueD, valueP, valueQ)

                self.__similarityMatrixD[currentPosition] = value

                # fill traceback matrix
                tracesP = []
                if valueP == case1P:
                    tracesP.append("upD")
                elif valueP == case2P:
                    tracesP.append("upP")

                tracesQ = []
                if valueQ == case1Q:
                    tracesQ.append("leftD")
                elif valueQ == case2Q:
                    tracesQ.append("leftQ")

                tracesD = []
                if value == valueD:
                    tracesD.append("diagD")
                if value == valueP:
                    tracesD.append("fromP")
                if value == valueQ:
                    tracesD.append("fromQ")

                self.__tracebackMatrixD[currentPosition] = tracesD
                self.__tracebackMatrixP[currentPosition] = tracesP
                self.__tracebackMatrixQ[currentPosition] = tracesQ

        return

    def traceback(self, complete_traceback):
        if complete_traceback:
            pass
        else:
            findingAlignment = True
            currentPosition = [self.__lengthA + 1, self.__lengthB + 1]
            traces = []
            while findingAlignment:
                if (currentPosition[0] == 0 or currentPosition[1] == 0):
                    findingAlignment = False
                    break
                else:
                    traces = self.__tracebackMatrixD[currentPosition]
                    traces = random.sample(traces, len(traces))
                    if traces[0] == "upD":
                        self.__alignA = self.__segA[currentPosition[0]]
                        self.__alignB = "-"
                        self.__alignment = " "
                        currentPosition[0] -= 1
                    elif traces[0] == "leftD":
                        self.__alignA = "-"
                        self.__alignB = self.__segB[currentPosition[1]]
                        self.__alignment = " "
                        currentPosition[1] -= 1

                    else:
                        self.__alignA = self.__segA[currentPosition[0]]
                        self.__alignB = self.__segB[currentPosition[1]]
                        if self.__alignA == self.__alignB:
                            self.__alignment = "*"
                        else:
                            self.__alignment = ":"
                        currentPosition[0] -= 1
                        currentPosition[1] -= 1
    def run(self,
            seq1_fasta_fn,
            seq2_fasta_fn,
            subst_matrix_fn,
            affine_cost_gap_open,
            affine_cost_gap_extend,
            complete_traceback):
            """Document me!"""
            self.__segA = seq1_fasta_fn
            self.__segB = seq2_fasta_fn
            self.__lengthA = len(seq1_fasta_fn)
            self.__lengthB = len(seq2_fasta_fn)
            self.__similarityMatrix = subst_matrix_fn
            self.__gapCost = affine_cost_gap_open
            self.__complete = complete_traceback

            self.initialize()
            self.fillMatrix()
            self.traceback(False)
            score = max(self.__similarityMatrix.values())
            return (score, self.__alignA, self.__alignment, self.__alignB)
