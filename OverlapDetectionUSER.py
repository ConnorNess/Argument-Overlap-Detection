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

################################################################################################################################
#If something is detected by all algorithms, might be worth noting what that is

def multipleOverlaps(multi_overlaps, levoverlaps, hamoverlaps, jarooverlaps, winkleroverlaps, jaccardoverlaps):
    
    dup_multi = []

    editoverlaps = itertools.product(levoverlaps, hamoverlaps)
    #Levenshtein and Hamming
    for a, b in editoverlaps:
        if a == b:
            dup_multi.append(a)

    ratiooverlaps = itertools.product(winkleroverlaps, jaccardoverlaps)
    #Winkler, and Jaccard
    for a, b in ratiooverlaps:
        if a == b:
            dup_multi.append(a)

    edit_and_ratio = [] #Everything but jaro

    #Get rid of duplicates
    for each in dup_multi:
        if each not in edit_and_ratio:
            edit_and_ratio.append(each)

    final_check = itertools.product(jarooverlaps, edit_and_ratio)
    #Finally, Jaro as our safe guard
    for a, b in final_check:
        if a == b:
            multi_overlaps.append(a)

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
    jaro_ratio = (match/ len(token1) + match / len(token2) + (match - transpositions + 1) / match)/ 3.0

    return(jaro_ratio)

def jaro_winkler(s1, s2, jaro_ratio):
    prefix = 0

    for i in range(min(len(s1), len(s2))):
        if (s1[i] == s2[i]) :
            prefix += 1
        else :
            break

    prefix = min(4, prefix) #Maximum 4 characters as indicated by the algorithms description
    winkler = jaro_ratio
    winkler += 0.1 * prefix * (1 - winkler)
    
    return(winkler)

def jaccard(token1, token2):
    token1chars = []
    for char in token1:
        token1chars.append(char)
    
    token2chars = []
    for char in token2:
        token2chars.append(char)

    intersection = len(list(set(token1chars).intersection(token2chars))) #AnB
    union = (len(set(token1chars)) + len(set(token2chars))) - intersection #AuB
    jaccard_similarity = (float(intersection) / union) #AnB / AuB
    return(jaccard_similarity)

################################################################################################################################
#Build array of arguments

def getArguments(arguments):
    #Read in the atom nodes of the supplied json files
    #Turn each json file into an array of strings holding the atoms + their ID
    
    domainpath = (os.path.dirname(os.path.realpath(__file__)) + '\\Pets\\SADFace\\Hand done\\') #<-------- Change directory here
    args = os.listdir(domainpath)
    

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
#Run

#Build arrays of atoms for each JSON
arguments = []
getArguments(arguments)

#Arrays for each overlap for each algorithm
levoverlaps = [] #All detected overlaps for an algorithm
hamoverlaps = []
jarooverlaps = []
jarowinkleroverlaps = []
jaccardoverlaps = []

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
    thisjaro_winkler = []
    #jaro + jaro winkler
    jaroratio = jaro(a_clean, b_clean)
    jaro_winkler_ratio = jaro_winkler(a_clean, b_clean, jaroratio)
    if (jaroratio >= 0.71): #ACCEPTABLE OVERLAP HERE
        thisjaro.append(a[1])
        thisjaro.append(b[1])
        jarooverlaps.append(thisjaro)
    if (jaro_winkler_ratio >= 0.8): #ACCEPTABLE OVERLAP HERE
        thisjaro_winkler.append(a[1])
        thisjaro_winkler.append(b[1])
        jarowinkleroverlaps.append(thisjaro_winkler)

    thisjaccard = []
    #jaccard
    jaccard_index = jaccard(a_clean, b_clean)
    if (jaccard_index >= 0.8): #ACCEPTABLE OVERLAP HERE
        thisjaccard.append(a[1])
        thisjaccard.append(b[1])
        jaccardoverlaps.append(thisjaccard)

multi_overlaps = []
multipleOverlaps(multi_overlaps, levoverlaps, hamoverlaps, jarooverlaps, jarowinkleroverlaps, jaccardoverlaps)

build_json(multi_overlaps, arguments)
print("done")