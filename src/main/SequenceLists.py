import sys
import fileinput

def SequenceLists(inputfile):
    sequences = {}
    currentSequence=""
    currentID =""
    content = []
    with open(inputfile) as f:
        content = f.readlines()
        for line in content:
            line = line.replace("\n","").replace("\r","")
            if line.startswith(">"):
                line.replace(">","")
                if currentSequence != "":
                    sequences[currentID] = currentSequence
                currentID = line
                currentSequence = ""
            elif line.isalpha():
                currentSequence = currentSequence + line

        sequences[currentID] = currentSequence
        f.close()

    return sequences