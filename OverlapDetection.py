import os
import re
import json
import itertools
from math import floor
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
    falsematches = 0
    matches = 0

    for links in handoverlaps: #each entry into the hand overlaps list
        for thisoverlap in overlaps: #each set of overlaps detected
            thismatch = 0 #tracks if both detected overlap nodes are present
            for handnodes in links: #each individual node in a list
                if thisoverlap[0] in handnodes or thisoverlap[1] in handnodes: 
                    thismatch += 1

            if (thismatch == 2):
                matches += 1

            if (thismatch == 1):
                falsematches += 1

    print("matches :")
    print(matches)
    print("false matches : ")
    print(falsematches)
                

################################################################################################################################
#ALGORITHMS HERE

def levenshtein(token1, token2):
    distances = np.zeros((len(token1) + 1, len(token2) + 1), dtype = int)

    for t1char in range(len(token1) + 1):
        distances[t1char][0] = t1char

    for t2char in range(len(token2) + 1):
        distances[0][t2char] = t2char
        
    a = 0
    b = 0
    c = 0
    
    for t1char in range(1, len(token1) + 1):
        for t2char in range(1, len(token2) + 1):
            if (token1[t1char-1] == token2[t2char-1]):
                distances[t1char][t2char] = distances[t1char - 1][t2char - 1]
            else:
                a = distances[t1char][t2char - 1]
                b = distances[t1char - 1][t2char]
                c = distances[t1char - 1][t2char - 1]
                
                if (a <= b and a <= c):
                    distances[t1char][t2char] = a + 1
                elif (b <= a and b <= c):
                    distances[t1char][t2char] = b + 1
                else:
                    distances[t1char][t2char] = c + 1
    
    return (distances[len(token1)][len(token2)])

def hamming(token1, token2):
    distance = 0
    for char in range (len(token2)): #loop for the shortest word
        if token1[char] != token2[char]:
            distance += 1
    distance += (len(token1) - len(token2)) #Add the remaining length of the input
    return distance

def jaro(token1, token2):
    maxdist = floor(max(len(token1), len(token2)) /2) -1
    match = 0
    hash1 = [0] * len(token1)
    hash2 = [0] * len(token2)

    for t1char in range(len(token1)):
        for t2char in range (max(0, t1char-maxdist), min(len(token2), t1char + maxdist+1)):

            if(token1[t1char] == token2[t2char] and hash2[t2char] == 0):
                hash1[t1char] = 1
                hash2[t2char] = 1
                match += 1
                break
    
    if(match == 0):
        return 0
    transpositions = 0
    point = 0

    for t1char in range (len(token1)):
        if(hash1[t1char]):
            while(hash2[point] == 0):
                point += 1
            if(token1[t1char] != token2[point]):
                point += 1
                transpositions += 1
    transpositions = floor(transpositions/2)

    return( (match / len(token1)) + (match / len(token2)) + ((match - transpositions+1)/match) /3.0) 

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
hamoverlaps = []
jarooverlaps = []

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

    thisham = []
    #hamming
    if len(a) < len(b):
        hamdist = hamming(b, a)
    else:
        hamdist = hamming(a, b)
    if(hamdist >= 5): #Again, the vibes dictate
        thisham.append(a[1])
        thisham.append(b[1])
        hamoverlaps.append(thisham)

    thisjaro = []
    #jaro
    jarodist = jaro(a[0], b[0])
    if (jarodist > 0.95): #Viiiiiiibes
        thisjaro.append(a[1])
        thisjaro.append(b[1])
        jarooverlaps.append(thisjaro)


# print("lev overlaps")
# print (levoverlaps)

# print("ham overlaps")
# print(hamoverlaps)

# print("jaro overlaps")
# print(jarooverlaps)

#for each algorithm array
print("levenshtein")
compareOverlaps(handoverlaps, levoverlaps)
print("____________________________")
print("hamming")
compareOverlaps(handoverlaps, hamoverlaps)
print("____________________________")
print("jaro")
compareOverlaps(handoverlaps, jarooverlaps)
print("____________________________")
#Total number of overlaps = every node in txt minus how many lines