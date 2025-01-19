close all;clear all;
clc
%% The program set1 is defined as the internal domain and set2 is defined as the SBPML domain
% 120 lines match the number of layers and 152 lines of material need to be edited in advance
% load *.txt at 21、22、23、26、32、35、82、91、105 lines
if exist('output1','dir')~=0
    rmdir('output1', 's')
end
if exist('output2','dir')~=0
    rmdir('output2', 's')
end
if exist('output3','dir')~=0
    rmdir('output3', 's')
end
pause(5)
dir_in1 = 'input';
dir_out1 = 'output1';
mkdir(dir_out1)

%%
node_all = dlmread([dir_in1,'\node_all.txt']);% abaqus model node
element_all = dlmread([dir_in1,'\element_all.txt']);% abaqus model elements
node_set1= importdata([dir_in1,'\node_set1.txt']); % set1 is the CPE4 elements
node_set1= node_set1(find( ~ isnan(node_set1))); 
node_set1 = sort(node_set1);
node_set2= importdata([dir_in1,'\node_set2.txt']); % Nodes of SBPML
node_set2= node_set2(find( ~ isnan(node_set2))); 
node_set2 = sort(node_set2);

% element_set1 = 61:220;%
% element_set1 = element_set1';
element_set1= importdata([dir_in1,'\element_set1.txt']); % CPE4 elements
element_set1= element_set1(find( ~ isnan(element_set1))); 
element_set1 = sort(element_set1);
element_set2= importdata([dir_in1,'\element_set2.txt']); % SBPML elements
element_set2= element_set2(find( ~ isnan(element_set2))); 
element_set2 = sort(element_set2);



for i = 1:length(node_set1)% Assign the corresponding coordinates to the abaqus node of the CPE4 elements
    m = find(node_all(:,1) == node_set1(i,1));
    node_set1(i,2:3) = node_all(m,2:3);
end

for i = 1:length(node_set2)% Assign coordinates to the abaqus node of the SBPML
    m = find(node_all(:,1) == node_set2(i,1));
    node_set2(i,2:3) = node_all(m,2:3);
end

for i = 1:length(element_set1)% abaqus connectivity for CPE4 elements
    m = find(element_all(:,1) == element_set1(i,1));
    element_set1(i,2:5) = element_all(m,2:5);
end
writematrix(element_set1,[dir_out1,'\element_set1.txt'])% CPE4 elements
for i = 1:length(element_set2)% abaqus connectivity for SBPML
    m = find(element_all(:,1) == element_set2(i,1));
    element_set2(i,2:5) = element_all(m,2:5);
end
% writematrix(element_set2,[dir_out1,'\element_set2.txt'])%abaqus SBPML elements
element_set1_fix = element_set1;
for i = 1:length(element_set1)% Replace CPE4 element abaqus connectivity with matlab connectivity
    for j = 2:5
        m2 = find(node_set1(:,1) == element_set1(i,j));
        element_set1_fix(i,j) = m2;
    end
end

element_set2_fix = element_set2;
for i = 1:length(element_set2)% Replace the abaqus connectivity for SBPML with matlab connectivity
    for j = 2:5
        m2 = find(node_set2(:,1) == element_set2(i,j));
        element_set2_fix(i,j) = m2;
    end
end
% txt output
writematrix(node_set1,[dir_out1,'\node.txt'])% matlab interior domain nodes
writematrix(node_set2,[dir_out1,'\node_pml.txt'])%matlab SBPML domain nodes
writematrix(element_set1_fix,[dir_out1,'\element.txt'])% matlab interior domain elements
writematrix(element_set2_fix,[dir_out1,'\element_pml.txt'])% matlab SBPML elements

% load_node= importdata([dir_in1,'\load_node.txt']); 
% load_node= load_node(find( ~ isnan(load_node))); 
% load_node = sort(load_node);
% for i = 1:length(load_node)
%     m = find(node_set1(:,1) == load_node(i,1));
%     load_node_fix(i,1) = m;
% end
% writematrix(load_node_fix,[dir_out1,'\Force.txt'])

cross= importdata([dir_in1,'\cross.txt']); % Load the node number that the SBPML shares with the CPE4 elements
cross= cross(find( ~ isnan(cross))); 
cross = sort(cross);
for i = 1:length(cross)% Replace the shared node number with the matlab interior domain node number
    m = find(node_set1(:,1) == cross(i,1));
    cross_set1(i,1) = m;
end
writematrix(cross_set1,[dir_out1,'\finite_inter.txt'])% matlab interior domain and SBPML domain intersection nodes number
for i = 1:length(cross)% Replace the shared node number with the matlab interior domain node number
    m = find(node_set2(:,1) == cross(i,1));
    cross_set2(i,1) = m;
end
writematrix(cross_set2,[dir_out1,'\infinite_inter.txt'])% matlab SBPML domain and interior domain intersection nodes number

fixnode= importdata([dir_in1,'\fixnode.txt']); % Fixed nodes number loaded in abaqus
fixnode= fixnode(find( ~ isnan(fixnode))); 
fixnode = sort(fixnode);
for i = 1:length(fixnode)% Replace the fixed nodes number in abaqus with the SBPML domain nodeS number
    m = find(node_set2(:,1) == fixnode(i,1));
    fixnode_fix(i,1) = m;
end
writematrix(fixnode_fix,[dir_out1,'\Artificial_boundary.txt'])% matlab fixed nodeS number
NNr = length(node_set1);
NEr = length(element_set1_fix);
NNpm = length(node_set2);
NEpm = length(element_set2_fix);
NNALL = length(node_all);
NEALL = length(element_all);

rows = 9;% The number of the Scaled Boundary Perfectly Matched Layers
ret1 = Only__ABC(NNr,NEr,NNpm,NEpm,rows);

[ret2,Inf_NOC_Nele2] = main2(NNALL,NEALL);

element_set2_2 = zeros(NEpm,5);
element_set2_2(:,1) = element_set2(:,1);
element_set2_2(:,2:5) = Inf_NOC_Nele2;

writematrix(element_set2_2,[dir_out1,'\element_set2.txt'])%abaqus SBPML elements

% m1= importdata('input/m1.txt');
% m1= m1(find( ~ isnan(m1)));
% save output3\m1.dat -ascii m1%
% 
% m2= importdata('input/m2.txt');
% m2= m2(find( ~ isnan(m2)));
% save output3\m2.dat -ascii m2%
% 
% m3= importdata('input/m3.txt');
% m3= m3(find( ~ isnan(m3)));
% save output3\m3.dat -ascii m3%

% for i = 1:7
%     eval(['m',char(string(i)),'=importdata("input\m',char(string(i)),'.txt");']);
%     eval(['m',char(string(i)),'= m',char(string(i)),'(find( ~ isnan(m',char(string(i)),')));']);
%     eval(['save output3\m',char(string(i)),'.dat -ascii m',char(string(i))]);
% end

%% 材料属性
IM1 = length(element_all);% Number of elements
IM2 = length(node_all);% Number of nodes
nmaterial = 2;% Number of types of material properties
for i = 1:nmaterial
    eval(['PML',char(string(i)),'=importdata("input\PML',char(string(i)),'.txt");']);
    eval(['PML',char(string(i)),'= PML',char(string(i)),'(find( ~ isnan(PML',char(string(i)),')));']);
    %eval(['save output3\m',char(string(i)),'.dat -ascii m',char(string(i))]);
end
materials = zeros(IM1,1);
for i = 1:nmaterial
    eval(['length1 = length(PML',char(string(i)),');']);
    for j = 1:length1
        eval(['materials(PML',char(string(i)),'(j,1),1) = i;']);
    end
end
save output3\materials.dat -ascii materials




