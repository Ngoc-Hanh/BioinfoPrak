import sys
import fileinput

def LoadSubstitutionMatix(matrixfile):
    matrix = {}
    key = ()
    value = 0
    path = sys.argv[0]
    dataDirectory = path.find("src")
    directory = path[0:dataDirectory-1] + "/data/"

    filename = ""
    if matrixfile == "pam":
        # load Pam matrix
        filename = directory + "pam250.txt"

    else:
        # load blosum matrix
        filename = directory + "blosum62.txt"

    content = []
    secondKeys = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M','F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
    with open(filename) as f:
        content = f.readlines()
        for line in content:
            if not (line.startswith("#") or (line.startswith(" ")) or (line.startswith("\n"))):
                #print(line)
                i = 0
                newLine = line.replace("  "," ")
                keys = newLine.split(" ")
                currentkey = keys[0]
                keys.remove(currentkey)
                for s in keys:
                    key = (currentkey, secondKeys[i])
                    matrix[key] = keys[i]
                    i += 1
                    #print(key, " - ", matrix [key])
            elif line.startswith(" "):
                #secondKeys = line.split(" ")
                #print(secondKeys[3])
                pass

        f.close()

    return matrix