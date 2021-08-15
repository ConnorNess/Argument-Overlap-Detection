import os
import re
import json
import itertools
from math import floor
import datetime
import numpy as np
import string
import uuid
import time

#Read the hand detected overlaps file by line creating arrays of each line that contributes to an overlap
#Read in the atom nodes of the supplied SADFace json files
#Run atoms through each algorithm comparing them to each other
#Flag atoms which are "close enough" as overlap
#Return these
#Check against hand detected overlaps
#Return these


################################################################################################################################
#output data to textfile

def build_json(overlaps, arguments):
    jsonoverlaps = []

    for thisoverlap in overlaps:

        overlap_id = str(uuid.uuid4())
        overlaps = []
        for thisnode in thisoverlap:
            for each in arguments:
                if thisnode == each[1]:
                    overlaps.append(each)
        
        overlap = {
            "overlap_id": overlap_id,
            "node_id_1": overlaps[0][0],
            "node_text_1": overlaps[0][1],
            "node_id_2": overlaps[1][0],
            "node_text_2": overlaps[1][1]
        }
        jsonoverlaps.append(overlap)

    filepath = os.path.dirname(os.path.realpath(__file__)) + '\\output'
    if not os.path.exists(filepath):
        os.makedirs(filepath)
    currenttime = datetime.datetime.now()
    filename = str(currenttime.day) + str(currenttime.hour) + str(currenttime.minute) + str(currenttime.second) + '.json'

    with open(os.path.join(filepath, filename), 'w') as f:
        json.dump(jsonoverlaps, f)

def output(multi_overlaps, multitrue, multitruecount, multifalse, multifalsecount, levoverlaps, levtrue, levtruecount, levfalse, levfalsecount, hamoverlaps, hamtrue, hamtruecount, hamfalse, hamfalsecount, jarooverlaps, jarotrue, jarotruecount, jarofalse, jarofalsecount):
    
    currenttime = datetime.datetime.now()
    filename = str(currenttime.day) + str(currenttime.hour) + str(currenttime.minute) + str(currenttime.second) + '.txt'

    filepath = os.path.dirname(os.path.realpath(__file__)) + '\\output'
    if not os.path.exists(filepath):
        os.makedirs(filepath)

    line_overlaps = 'Detected by Levenshtein, Hamming, and Jaro:'
    with open(os.path.join(filepath, filename), 'w') as f:
        f.write('Detected by Levenshtein, Hamming, and Jaro:\n')
        f.write('\nCorrect Detection: ' + str(multitruecount) + '\n')
        for overlap in multitrue:
            f.write(str(overlap) + '\n')
        f.write('\nIncorrect Detection: ' + str(multifalsecount) + '\n')
        for overlap in multifalse:
            f.write(str(overlap) + '\n')
        f.write('\n____________________\n\n')

        f.write('\nDetected by Levenshtein:\n')
        f.write('\nCorrect Detection: ' + str(levtruecount) + '\n')
        for overlap in levtrue:
            f.write(str(overlap) + '\n')
        f.write('\nIncorrect Detection: ' + str(levfalsecount) + '\n')
        for overlap in levfalse:
            f.write(str(overlap) + '\n')
        f.write('\n____________________\n\n')

        f.write('\nDetected by Hamming:\n')
        f.write('\nCorrect Detection: ' + str(hamtruecount) + '\n')
        for overlap in hamtrue:
            f.write(str(overlap) + '\n')
        f.write('\nIncorrect Detection: ' + str(hamfalsecount) + '\n')
        for overlap in hamfalse:
            f.write(str(overlap) + '\n')
        f.write('\n____________________\n\n')

        f.write('\nDetected by Jaro:\n')
        f.write('\nCorrect Detection: ' + str(jarotruecount) + '\n')
        for overlap in jarotrue:
            f.write(str(overlap) + '\n')
        f.write('\nIncorrect Detection: ' + str(jarofalsecount) + '\n')
        for overlap in jarofalse:
            f.write(str(overlap) + '\n')
        f.write('\n____________________\n\n')

################################################################################################################################
#If something is detected by all algorithms, might be worth noting what that is

def multipleOverlaps(multi_overlaps, overlapsA, overlapsB, overlapsC):
    
    overlaps = itertools.product(overlapsA, overlapsB, overlapsC)

    for a, b, c in overlaps:
        if a == b == c:
            multi_overlaps.append(a)

################################################################################################################################
#Compare detected overlaps to self detected overlaps

def compareOverlaps(handoverlaps, overlaps, correctoverlap, correctcount, incorrectoverlap, incorrectcount):
    

    for overlap in overlaps:
        for handoverlap in handoverlaps:

            match = 0

            if overlap[0] == handoverlap[0] or overlap[0] == handoverlap[1]:
                match += 1
            if overlap[1] == handoverlap[0] or overlap[1] == handoverlap[1]:
                match += 1
            if overlap[0] == overlap[1]: #this *should* never hit thanks to itertools but, better safe
                match = 0

            if match == 2:
                break

        
        if match == 2:
            correctoverlap.append(overlap)
            correctcount += 1
        else:
            incorrectoverlap.append(overlap)
            incorrectcount += 1

    return correctcount, incorrectcount


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
    chars_t1 = [0] * len(token1)
    chars_t2 = [0] * len(token2)

    if (token1 == token2):
        return 1

    for t1char in range(len(token1)):
        for t2char in range (max(0, t1char-maxdist), min(len(token2), t1char + maxdist+1)):

            if(token1[t1char] == token2[t2char] and chars_t2[t2char] == 0):
                chars_t1[t1char] = 1
                chars_t2[t2char] = 1
                match += 1
                break
    
    if(match == 0):
        return 0
    transpositions = 0
    point = 0

    for t1char in range (len(token1)):
        if(chars_t1[t1char]):
            while(chars_t2[point] == 0):
                point += 1
            if(token1[t1char] != token2[point]):
                point += 1
                transpositions += 1
    transpositions = floor(transpositions/2)
    jaroratio = (match/ len(token1) + match / len(token2) + (match - transpositions + 1) / match)/ 3.0

    return(jaroratio) 

def synonyms_antonyms(token1, token2):

    

    return(0)

################################################################################################################################
#Build array of arguments

def getArguments(arguments):
    #Read in the atom nodes of the supplied AIF json files
    #Turn each json file into an array of strings holding the atoms + their ID
    
    domainpath = (os.path.dirname(os.path.realpath(__file__)) + '\\Pets\\SADFace\\Hand done\\') #<-------- Change directory here
    #domainpath = (os.path.dirname(os.path.realpath(__file__)) + '\\welfare\\this\\')
    args = os.listdir(domainpath)
    #args = ['2229.json', '2235.json', 'self_1.json'] #Used for testing quickly
    

    for each in args:
            with open(domainpath + each) as arg:

                data = json.load(arg) #JSON to dict
                nodes = data.get('nodes') #Dict to list
                count = len(nodes) #How many nodes long is this json?

                for i in range(count):
                    thisarg = [] #temp array to hold text and id, will append once filled
                    try:
                        thisarg.append(data['nodes'][i]['text']) #Both SADFace and AIFdb's JSON output of AIF use 'text' under 'nodes', makes this a bit easier
                        try:
                            thisarg.append(data['nodes'][i]['id'])
                            arguments.append(thisarg) #Only append an argument with both text and id
                        except KeyError:
                            try:
                                if(data['nodes'][i]['type'] == "I"): #if this is AIF, we only want i-nodes
                                    thisarg.append(data['nodes'][i]['nodeID']) 
                                    arguments.append(thisarg) #Only append an argument with both text and id
                            except KeyError:
                                print("no id") #uh oh, somethin wrong with the json
                    except KeyError: #Means no text is present - its a scheme node likely
                        pass #If it ain't got text, we got no use for it

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
getArguments(arguments) #<-----------Change type here (AIF or SADFace)

#Arrays for each overlap for each algorithm
levoverlaps = [] #All detected overlaps for an algorithm
hamoverlaps = []
jarooverlaps = []


#For each atom in a SADFACE
#Compare to each atom in all SADFaces
#Using itertools combinations to easily compare all elements but only once, two for loops would double up
for a, b in itertools.combinations(arguments, 2):

    a_no_punctuation = a[0].translate(str.maketrans('', '', string.punctuation))
    a_clean = a_no_punctuation.lower()
    b_no_punctuation = b[0].translate(str.maketrans('', '', string.punctuation))
    b_clean = b_no_punctuation.lower()

    thislev = []
    #levenshtein
    levdist = levenshtein(a_clean, b_clean) #0 = the text, 1 = the id
    if (levdist <= 12): #ACCEPTABLE OVERLAP HERE - gotten through testing
        thislev.append(a[1])
        thislev.append(b[1])
        levoverlaps.append(thislev)

    thisham = []
    #hamming
    if len(a_clean) < len(b_clean):
        hamdist = hamming(b_clean, a_clean)
    else:
        hamdist = hamming(a_clean, b_clean)
    if(hamdist <= 13): #ACCEPTABLE OVERLAP HERE
        thisham.append(a[1])
        thisham.append(b[1])
        hamoverlaps.append(thisham)

    thisjaro = []
    #jaro
    jaroratio = jaro(a_clean, b_clean)
    if (jaroratio >= 0.71): #ACCEPTABLE OVERLAP HERE
        thisjaro.append(a[1])
        thisjaro.append(b[1])
        jarooverlaps.append(thisjaro)

levtrue = [] #For all correctly identified overlaps
levtruecount = 0
levfalse = [] #For all incorrecntly identified overlaps
levfalsecount = 0
levtruecount, levfalsecount = compareOverlaps(handoverlaps, levoverlaps, levtrue, levtruecount, levfalse, levfalsecount)

hamtrue = []
hamtruecount = 0
hamfalse = []
hamfalsecount = 0
hamtruecount, hamfalsecount = compareOverlaps(handoverlaps, hamoverlaps, hamtrue, hamtruecount, hamfalse, hamfalsecount)

jarotrue = []
jarotruecount = 0
jarofalse = []
jarofalsecount = 0
jarotruecount, jarofalsecount = compareOverlaps(handoverlaps, jarooverlaps, jarotrue, jarotruecount, jarofalse, jarofalsecount)


multi_overlaps = []
multipleOverlaps(multi_overlaps, levoverlaps, hamoverlaps, jarooverlaps)

multitrue = []
multitruecount = 0
multifalse = []
multifalsecount = 0
multitruecount, multifalsecount = compareOverlaps(handoverlaps, multi_overlaps, multitrue, multitruecount, multifalse, multifalsecount)


#Oops, all meaningful data! - this looks messy~~~~
output(multi_overlaps, multitrue, multitruecount, multifalse, multifalsecount, levoverlaps, levtrue, levtruecount, levfalse, levfalsecount, hamoverlaps, hamtrue, hamtruecount, hamfalse, hamfalsecount, jarooverlaps, jarotrue, jarotruecount, jarofalse, jarofalsecount)
build_json(multi_overlaps, arguments)
print("done")