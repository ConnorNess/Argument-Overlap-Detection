# Argument-Overlap-Detection
Detection of overlapping argument nodes in SADFace

FOR GENERAL USE </br>
Replace the domainpath variable within def getArguments to a directory consisting ONLY of JSON files. </br>
All atom nodes and their IDs will be retireved and each will be compared to all other nodes through several algorithms </br>
These algorithms will output their detected overlaps by ID </br>
 </br>
FOR TESTING/DEVELOPMENT </br>
Self detected overlaps are listed within a txt file to fill lists, each line within the txt file are the overlapping nodes. </br>
The program will be run the same as general use. </br>
After these algorithms have filled a list output, the lists are compared to hand detected overlaps with def compareOverlaps. </br>


This work was done for the purpose of a dissertation for a Masters of Science in Computing @ Edinburgh Napier. Interested in what I do or are hiring? Send a message and I'll reply <3
