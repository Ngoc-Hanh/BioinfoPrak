from InOutManager.OutputBase import OutputBaseClass
import math

@OutputBaseClass.register
class PairwiseOutputManager(OutputBaseClass):
    """
    Base class for displaying or saving output files in correct format.
    """
    def aligned_character(self, align1, align2):
        alignment = ""
        for i in range(len(align1)):
            if align1[i] == align2[i]:
                alignment += "*"
            elif align1[i] == "-" or align2[i] == "-":
                alignment += "0"
            else:
                alignment += ":"

        return alignment

    def display(self, score, result1, result2):
        '''
        display output
        '''
        #id1, id2, score, aligned_sequence_list
        # alignment_result = list(*result)
        # score = alignment_result[2]
        # res= list(alignment_result[3])
        # result1 = res[0]
        # result2 = res[1]
        alignment = self.aligned_character(result1, result2)
        align = "".join(c for c in alignment)
        line_counter = math.floor(len(result1)/80)
        print("optimal alignment score : {}".format(score))
        counter = 0
        while counter <= line_counter:
            print(result1[counter*line_counter:counter*line_counter + 80])
            print(str(alignment[counter * line_counter:counter * line_counter + 80]))
            print(result2[counter*line_counter:counter*line_counter + 80])
            counter += 1

    def save(self, score, result1, result2, filepath):
        '''
        save output
        :param filepath: path to the output file
        '''
        # alignment_result = tuple(*result)
        # score = alignment_result(2)
        # result1= alignment_result(3)
        # result2 = alignment_result(4)
        alignment = self.aligned_character(result1, result2)
        line_counter = math.floor(len(result1)/80)
        counter = 0
        filename = filepath + "\\" + "result.fas"
        with open(filename, 'w') as f:
            print("optimal alignment score : {}".format(score))
            while counter <= line_counter:
                f.write(result1[counter*line_counter:counter*line_counter + 80])
                f.write(alignment[counter * line_counter:counter * line_counter + 80])
                f.write(result2[counter*line_counter:counter*line_counter + 80])
                counter +=1
        f.close()

@OutputBaseClass.register
class FastaOutputManager(OutputBaseClass):

    """
    Base class for displaying or saving output files in correct format.
    """

    def display(self, score, id1, result1, id2, result2, id3, result3):
        '''
        Validate the content of the input file. Exit with error if the file is corrupted or contains strange character
        :param filepath: path to the input file
        :return: TRUE if the file can be loaded with the correct format
        '''
        #id1, id2, score, aligned_sequence_list
        # alignment_result = list(*result)
        # score = alignment_result[2]
        # res= list(alignment_result[3])
        # result1 = res[0]
        # result2 = res[1]

        line_counter = math.floor(len(result1)/80)
        print("optimal alignment score : {}".format(score))
        counter = 0
        while counter <= line_counter:
            print (id1, end="  ")
            print(result1[counter*line_counter:counter*line_counter + 80])
            print (id2, end="  ")
            print(result2[counter*line_counter:counter*line_counter + 80])
            print (id3, end="  ")
            print(result3[counter*line_counter:counter*line_counter + 80])

            counter += 1

    def save(self, score, id1, result1, id2, result2, id3, result3, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the output file
        '''
        # alignment_result = tuple(*result)
        # score = alignment_result(2)
        # result1= alignment_result(3)
        # result2 = alignment_result(4)
        alignment = self.aligned_character(result1, result2)
        line_counter = math.floor(len(result1)/80)
        counter = 0
        filename = filepath + "\\" + "result.fas"
        with open(filename, 'w') as f:
            print("optimal alignment score : {}".format(score))
            while counter <= line_counter:
                f.write(id1, end="  ")
                f.write(result1[counter * line_counter:counter * line_counter + 80])
                f.write(id2, end="  ")
                f.write(result2[counter * line_counter:counter * line_counter + 80])
                f.write(id3, end="  ")
                f.write(result3[counter * line_counter:counter * line_counter + 80])

                counter +=1
        f.close()

@OutputBaseClass.register
class NewickOutputManager(OutputBaseClass):

    """
    Base class for displaying or saving output files in correct format.
    """

    def display(self, result):
        '''
        Validate the content of the input file. Exit with error if the file is corrupted or contains strange character
        :param filepath: path to the input file
        :return: TRUE if the file can be loaded with the correct format
        '''
        pass

    def save(self, result, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a dictionary of inputs
        '''
        pass

@OutputBaseClass.register
class RNAOutputManager(OutputBaseClass):

    """
    Base class for displaying or saving output files in correct format.
    """

    def display(self, result):
        '''
        Validate the content of the input file. Exit with error if the file is corrupted or contains strange character
        :param filepath: path to the input file
        :return: TRUE if the file can be loaded with the correct format
        '''
        pass

    def save(self, result, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a dictionary of inputs
        '''
        pass