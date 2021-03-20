# Structural Calculators

## How to use "Truss.m"
Open the file "Truss.m" input:

0. Decide whether you would like to calculate internal forces or not (internal forces could only be calculated on truss structures)
1. Input nodes coordinates in "nodes"
2. Input elements, from node a to node b in "elements"
3. Input support type, i.e. how many unkown reactions are there on each node, in "SupportTypesOnNodes"
4. Input external forces or moments in "ExF" and "ExM" respectively
5. Run
6. Obtained results: internal forces on each element, and reaction forces on each supporting nodes.

## How to use "truss_DSM.m"
Open the file "truss_DSM.m" input:
1. Input nodes coordinates in "nodes"
2. Input elements, from node a to node b in "ele"
3. Input geometry properties, which is a vector that contains properties of each element
4. Input support types, 2 means both x and y are constrained; 12 means y constrained, 11 means x constrained, in "Supporting"
5. Input external forces or moments in "ExF" and "ExM" respectively
6. Run
7. Obtained results: displacement of each degrees, internal forces on each element, and reaction forces on each supporting nodes.

## How to use "frame_DSM.m"
Open the file "frame_DSM.m" input:
1. Input nodes coordinates in "nodes"
2. Input elements, from node a to node b in "eles"
3. Input geometry properties, which is a vector that contains properties of each element
4. Input support types, as instructed by the comments
5. Input external nodal forces and element equivalent forecs
6. Run
7. Obtained results: displacement of each degrees, element end forces, and reaction forces on each supporting nodes, and a simple visualization of the structure, you can change the variable "mag_factor" to exaggerate the deformation.

“resolution” relates to the level of finess of visualization. If you would like to check the element end forces, change "resolution" to 0.

## How to use "Frame.m"
Similar with "Truss", but you can input load in this script. See annotation in the script.
Obtained results: reaction forces on each supporting nodes.

## Demostration on "Truss.m"
![image](demo_problem.png)
![image](demo_problem_result.jpg)
![image](demo_problem_result_figure.jpg)

And "truss_DSM.m" can provide exactly the same result (but no visualization).

## Demostration on "frame_DSM.m"
![image](demo_frame_config.png)
![image](demo_frame_result.png)
![image](demo_frame_fig.png)
