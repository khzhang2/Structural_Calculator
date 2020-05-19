# Structural Calculators

## how to use "Truss"
0. choose whether you would like to calculate internal forces or not (internal forces could only be calculated on truss structures)
1. input nodes coordinates in "nodes"
2. input elements, from node a to node b in "elements"
3. input support type, i.e. how many unkown reactions are there on each node, in "SupportTypesOnNodes"
4. input external forces or moments in "ExF" and "ExM" respectively
5. run
6. obtained results: internal forces on each element, and reaction forces on each supporting nodes.

## how to use "Frame"
Similar with "Truss", but you can input load in this script. See annotation in the script.
Obtained results: reaction forces on each supporting nodes.
