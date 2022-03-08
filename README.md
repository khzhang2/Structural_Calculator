# Structural Calculators
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

## Demostration on "Truss.m"
![image](demo_problem.png=100x)
![image](demo_problem_result.jpg=100x)
![image](demo_problem_result_figure.jpg=100x)

And "truss_DSM.m" can provide exactly the same result (but no visualization).

## Demostration on "frame_DSM.m"
![image](demo_frame_config.png)
![image](demo_frame_result.png)
![image](demo_frame_fig.png)
