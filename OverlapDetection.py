import os
import re
import json
import itertools
import numpy as np

#Read the hand detected overlaps file by line creating arrays of each line that contributes to an overlap
#Read in the atom nodes of the supplied SADFace json files
#Run atoms through each algorithm comparing them to each other
#Flag atoms which are "close enough" as overlap
#Return these
#Check against hand detected overlaps
#Return these


################################################################################################################################
#Compare detected overlaps to self detected overlaps

def compareOverlaps(handoverlaps, overlaps):
    #Check against hand detected overlaps

    for links in handoverlaps: #each entry into the hand overlaps list
        for thisoverlap in overlaps: #each set of overlaps detected
            matches = 0 #tracks if both detected overlap nodes are present
            for handnodes in links: #each individual node in a list
                if thisoverlap[0] in handnodes or thisoverlap[1] in handnodes: 
                    matches += 1

            if (matches == 2):
                print(thisoverlap)

################################################################################################################################
#ALGORITHMS HERE

def levenshtein(token1, token2):
    distances = np.zeros((len(token1) + 1, len(token2) + 1), dtype = int)

    for char1 in range(len(token1) + 1):
        distances[char1][0] = char1

    for char2 in range(len(token2) + 1):
        distances[0][char2] = char2
        
    a = 0
    b = 0
    c = 0
    
    for char1 in range(1, len(token1) + 1):
        for char2 in range(1, len(token2) + 1):
            if (token1[char1-1] == token2[char2-1]):
                distances[char1][char2] = distances[char1 - 1][char2 - 1]
            else:
                a = distances[char1][char2 - 1]
                b = distances[char1 - 1][char2]
                c = distances[char1 - 1][char2 - 1]
                
                if (a <= b and a <= c):
                    distances[char1][char2] = a + 1
                elif (b <= a and b <= c):
                    distances[char1][char2] = b + 1
                else:
                    distances[char1][char2] = c + 1
    
    return (distances[len(token1)][len(token2)])


################################################################################################################################
#Build array of arguments

def getSADFaces(arguments):
    #Read in the atom nodes of the supplied SADFace json files
    #Turn each json file into an array of strings holding the atoms + their ID

    domainpath = (os.path.dirname(os.path.realpath(__file__)) + "\\Pets\\SADFace\\Hand Done\\")
    SADFaces = os.listdir(domainpath)
    #SADFaces = ['2229.json', '2235.json'] #Used for testing quickly

    for each in SADFaces:
            with open(domainpath + each) as SADFace:

                data = json.load(SADFace) #JSON to dict
                nodes = data.get('nodes') #Dict to list
                count = len(nodes) #How many nodes long is this json?

                for i in range(count):
                    thisarg = [] #temp array to hold text and id, will append once filled
                    try:
                        thisarg.append(data['nodes'][i]['text'])
                        try:
                            thisarg.append(data['nodes'][i]['id'])
                            arguments.append(thisarg) #Only append an argument with both text and id
                        except KeyError:
                            print("no id") #uh oh, somethin wrong with the json
                    except KeyError: #Means no text is present - its a scheme node likely
                        pass

################################################################################################################################
#Build array of self detected overlaps

def getHandOverlaps(handoverlaps):
    #Read the hand detected overlaps file by line creating arrays of each line that contributes to an overlap

    handoverlapsfilepath = os.path.dirname(os.path.realpath(__file__)) + '\\Pets\\overlaps.txt'
    with open(handoverlapsfilepath) as handoverlapsfile:
        handoverlaps.extend([line.split() for line in handoverlapsfile])

    handoverlapsfile.close()

################################################################################################################################
#Run

#Build array using getHandOverlaps
handoverlaps = []
getHandOverlaps(handoverlaps)
#Build arrays of atoms for each SADFace
arguments = []
getSADFaces(arguments)

#Arrays for each overlap for each algorithm
levoverlaps = []

#For each atom in a SADFACE
#Compare to each atom in all SADFaces
#Using itertools combinations to easily compare all elements but only once, two for loops would double up
for a, b in itertools.combinations(arguments, 2):
    
    thislev = []
    #levenshtein
    levdist = levenshtein(a[0],b[0]) #0 = the text, 1 = the id
    if (levdist <= 5): #Ok, kinda have to just chose a random-ish number based on what the vibe is - testing will show whats better FOR THIS DATASET but something more algorithmic should replace just a flat integer value
        thislev.append(a[1])
        thislev.append(b[1])
        levoverlaps.append(thislev)

    #hamming
        #if match add to hamming array
    #etc.
        #if match add to x array

#for each algorithm array

compareOverlaps(handoverlaps, levoverlaps)

#Total number of overlaps = every node in txt minus how many lines