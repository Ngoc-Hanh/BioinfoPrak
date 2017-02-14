from InOutManager.InputBase import InputBaseClass

class InputManager:
    """
    Bases for validating and loading input files.
    """
    amino_acid =['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M','F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

    def __init__(self):
        self.not_valid_sequences = False
        self.not_valid_substitution_matrix = False

    def load_fasta_sequences(self, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a dictionary of inputs
        '''
        sequences = {}
        current_sequence = ""
        current_id = ""
        content = []
        with open(filepath) as f:
            content = f.readlines()
            for line in content:
                if line.startswith(";"): continue
                line = line.replace("\n", "").replace("\r", "")
                if line.startswith(">"):
                    line.replace(">", "")
                    if current_sequence != "":
                        sequences[current_id] = current_sequence
                    current_id = line
                    current_sequence = ""
                else:
                    for c in list(line):
                        if c not in self.amino_acid:
                            self.not_valid_sequences = True
                            break
                            return None
                    current_sequence = current_sequence + line

        sequences[current_id] = current_sequence
        f.close()

        return sequences

    def load_substitution_matrix(self, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a dictionary of inputs
        '''
        matrix = {}
        key = ()
        value = 0

        content = []
        secondKeys = []
        #['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
        with open(filepath) as f:
            content = f.readlines()
            for line in content:
                if not (line.startswith("#") or (line.startswith(" ")) or (line.startswith("\n"))):
                    # print(line)
                    i = 0
                    newLine = line.replace("  ", " ")
                    keys = newLine.split(" ")
                    currentkey = keys[0]
                    keys.remove(currentkey)
                    for s in keys:
                        key = (currentkey, secondKeys[i])
                        matrix[key] = keys[i]
                        i += 1
                        # print(key, " - ", matrix [key])
                elif line.startswith("   "):
                    secondKeys = line.replace("   ","").replace("  "," ").replace("\n","").split(" ")
                    for k in secondKeys:
                        if k not in self.amino_acid:
                            self.not_valid_substitution_matrix = True
                            break
                            return
                    #print(secondKeys[3])

            f.close()

        return matrix

class RNAInputSequence():
    """
    Base class for validating and loading input files.
    """

    def validate(self, filepath):
        '''
        Validate the content of the input file. Exit with error if the file is corrupted or contains strange character
        :param filepath: path to the input file
        :return: TRUE if the file can be loaded with the correct format
        '''
        pass

    def load(self, filepath):
        '''
        Load the input file for future uses
        :param filepath: path to the input file
        :return: a dictionary of inputs
        '''
        pass