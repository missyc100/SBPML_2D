------------------------------------------------------------Quick Guide to Running Our Cases---------------------------------------------------------------------------

1. Modify the file reading path in U1004-E24.for (refer to lines 13-17 in U1004-E24.for)
2. Create a job for Job-PML.inp file in ABAQUS and associate it with U1004-E24.for

-----------------------------The tutorial is divided into two parts: Part 1 is SBPML UEL element, Part 2 is Domain Reduction Method----------------------------------------------------------

The environment I use is ABAQUS 2022, MATLAB 2022a, Intel(R) MPI Library 2019 Update 9
The main program in the "Calculated wave velocity" folder is main.m
The main program in the "Generate matching layer parameters with one click" folder is main1.m
The main program in the DRM folder is main.m

1. Install the Abaqus2Matlab.mlappinstall plugin in MATLAB
2. Run Calculated wave velocity\main to get the recommended element size, matching layer thickness, and absorption coefficient
3. Create a model similar to model1.cae. The SBPML domain needs to be created in PART1, each layer of the SBPML domain must be kept parallel, select implicit dynamics for analysis step, and select asymmetric for matrix storage method
4. Create the following sets in PART1:
   - cross: Node set at the interface between SBPML elements and CPE4 elements
   - element_all: Set of all elements in PART1
   - element_set1: Set of CPE4 elements
   - element_set2: Set of SBPML elements
   - fixnode: Fixed node set
   - node_all: Set of all nodes in PART1
   - node_set1: Set of CPE4 nodes
   - node_set2: Set of SBPML nodes
   - PML1: Set of SBPML elements for material 1
   - PML2: Set of SBPML elements for material 2
   - m1: Set of CPE4 elements for material 1
   - m2: Set of CPE4 elements for material 2
5. Export the INP file from CAE and copy the corresponding named sets to '\Generate matching layer parameters with one click\input' path
6. Modify line 120 (matching layer number) and line 152 (number of materials) in \Generate matching layer parameters with one click\main1.m. In this tutorial, there are 2 materials
7. Run main1.m
8. Create a 'read' folder in the CAE file directory, copy all files from output3 (generated in step 7) to the read folder
9. Modify the INP file exported in step 5:
   (1) Delete element connectivity in PART1
   (2) Copy the UEL definition format to the appropriate position (refer to line 813 in Job-PML.inp)
   (3) Modify the material properties assigned to CPE4 (refer to lines 1092, 1095 in Job-PML.inp)
10. Modify the file reading path, number of element nodes, material properties, and absorption parameters in U1004-E24.for (refer to line 2 in U1004-E24.for)
11. Create a job for Job-PML.inp file in ABAQUS and associate it with U1004-E24.for

-----------------------------------------SBPML is now working properly, below is the seismic input section---------------------------------------------------------------------

1. From the INP file exported in step 5 above, copy the corresponding named sets to '\DRM\input' path
2. Select a PART for DRM input (in this example, PART1 is selected), export PART1 to a new INP file, and add commands to output mass and stiffness matrices at the end of the INP file (refer to lines 1641-1650 in Job-KKMM.inp)
3. Run the INP file from step 2 in ABAQUS. After completion, you will get two *.mtx files. Copy these files to '\DRM\MTX' path
4. Determine the desired seismic time history input (refer to line 63 in '\DRM\main.m')
5. Determine the type of wave input: options are 1-P wave, 2-SV wave (refer to line 72 in '\DRM\main.m'). Set the seismic incidence angle (refer to line 65 in '\DRM\main.m'). Note that Fei1 is adjustable, Fei2 is 0
6. Modify bedrock material properties (refer to lines 8-10 in '\DRM\OneDimension_FE_AB.m')
7. Modify stratum material properties (refer to '\DRM\material_constant_Solid.txt')
   Note: You need to check the content of matrix YY when '\DRM\main.m' runs to line 72, this will tell you how many layers of materials need to be defined
   Do not delete the lamd content in the first line of '\DRM\material_constant_Solid.txt'
8. Create a case1 folder in the CAE file directory, copy the contents of \DRM\output1, \DRM\output2, and the \DRM\output4 folder to the case1 folder (refer to \DRM\case1)
9. Include the txt files obtained in step 8 in the INP file (refer to lines 2404, 2446, 2476 in Job-PML-DRM.inp)
10. Create a job for Job-PML-DRM.inp file in ABAQUS and associate it with U1004-E24.for

---------------------------------------------------------End-----------------------------------------------------------------------------------------------------------