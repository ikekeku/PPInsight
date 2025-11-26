This folder contains sample output from the Lightdock simple protein-protein docking tutorial (https://lightdock.org/tutorials/0.9.3/simple_docking) for the purpose of buildling Task #3 of PPInsight.
'simulation/swarm_0' is a folder containg expected example output files from running the simulation in the tutorial linked above after correctly installing Rosetta (https://github.com/lightdock/lightdock/tree/master?tab=readme-ov-file#3-installation). Some folders/files have been omitted to make navigation simpler (i.e., swarm_(11-19 and 21-181):

**The values we want are in the 'gso_X.out' files under 'Scoring'** (where 'X' is the step number):
> "If the [simulation] run is successful... we will find the output files for each of the independent swarms. In this case since we only simulated 1 swarm, the results will be inside swarm_0. The output files will be named as gso_X.out, being X is the step number.
> **NOTE** LightDock will only store the results every 10 simulation steps (0, 10, 20, 30, etc.)
> In each of the output files, every line corresponds to a glowworm agent in the GSO algorithm. The numbers enclosed by ( and ), refer to the [x,y,z] coordinates in the translational space and the quaternion vector (q = a + 0i + 0j + 0k) in the rotational space. If ANM is enabled, this vector would expand by the number of normal modes considered for receptor and ligand respectively.
> The coordinates of the predicted complexes are followed by the ID of the complex and the last column refers to the scoring, in this case as calculated with DFIRE:

> ```bash head -2 swarm_0/gso_100.out ```

> ```#Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring ```
> 
> ``` ( 8.3619938, -9.0568297, -23.0839562, -0.6368250, -0.1214085,  0.4488531, -0.6150161,  5.2922640,  4.6293871,  5.1679589,  3.6238049,  3.9717360,  7.2390017,  4.5821031,  5.4562578,  1.9734713,  2.9626742,  1.0323748,  5.8207747,  3.4377706,  2.6552282,  4.8705878,  3.0931470,  3.4670304,  3.9487060,  2.4948334,  4.5626560)    0    0  27.83190014  6 0.320  18.65982945 ```

**Basics of LightDock (very helpful):** https://lightdock.org/tutorials/0.9.3/basics
