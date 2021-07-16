# SADFace-Overlap-Detection
Detection of overlapping argument nodes in SADFace

FOR GENERAL USE
Replace the domainpath variable within def getSADFaces to a directory consisting ONLY of SADFace JSON files.
All atom nodes and their IDs will be retireved and each will be compared to all other nodes through several algorithms
These algorithms will output their detected overlaps by ID

FOR TESTING/DEVELOPMENT
Self detected overlaps are listed within a txt file to fill lists, each line within the txt file are the overlapping nodes.
The program will be run the same as general use.
After these algorithms have filled a list output, the lists are compared to hand detected overlaps with def compareOverlaps
