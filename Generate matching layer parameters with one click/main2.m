function [ret2,Inf_NOC_Nele2] = main2(NNALL,NEALL)
% close all;clear all;
% clc
ele_first = dlmread('output2\ele_first.dat');
Inf_NOC_Na = dlmread('output2\Inf_NOC_Na.dat');
Inf_NOC_Nb = dlmread('output2\Inf_NOC_Nb.dat');
Inf_NOC_Nele = dlmread('output2\Inf_NOC_Nele.dat');
NOC_pm = dlmread('output2\NOC_pm.dat');
XofN_pm = dlmread('output2\XofN_pm.dat');
element_set2_fix = dlmread('output1\element_pml.txt');
node_set2_fix = dlmread('output1\node_pml.txt');

node_all = dlmread('input\node_all.txt',',',[0,0,NNALL-1,2]);
element_all = dlmread('input\element_all.txt',',',[0,0,NEALL-1,4]);

m = length(element_all);%216
n = length(node_all);%259
ele_first_fix = zeros(m,2);
Inf_NOC_Na_fix = zeros(m,2);
Inf_NOC_Nb_fix = zeros(m,2);
Inf_NOC_Nele_fix = zeros(m,4);
NOC_pm_fix = zeros(m,4);
XofN_pm_fix = zeros(n,2);

for i = 1:length(ele_first)
    for j = 1:2
        Inf_NOC_Nb2(i,j) = node_set2_fix(Inf_NOC_Nb(i,j),1);
        Inf_NOC_Na2(i,j) = node_set2_fix(Inf_NOC_Na(i,j),1);
    end
    for k = 1:4
        Inf_NOC_Nele2(i,k) = node_set2_fix(Inf_NOC_Nele(i,k),1);
        NOC_pm2(i,k) = node_set2_fix(NOC_pm(i,k),1);
    end
end

for i = 1:length(ele_first)
    ele_first_fix(element_set2_fix(i,1),:) = ele_first(i,:);
    Inf_NOC_Na_fix(element_set2_fix(i,1),:) = Inf_NOC_Na2(i,:);
    Inf_NOC_Nb_fix(element_set2_fix(i,1),:) = Inf_NOC_Nb2(i,:);
    Inf_NOC_Nele_fix(element_set2_fix(i,1),:) = Inf_NOC_Nele2(i,:);
    NOC_pm_fix(element_set2_fix(i,1),:) = NOC_pm2(i,:);

end

for i = 1:length(XofN_pm)
    XofN_pm_fix(node_set2_fix(i,1),:) = XofN_pm(i,:);
end
dir = 'output3';
mkdir(dir)

save output3\ele_first.dat -ascii ele_first_fix
save output3\Inf_NOC_Na.dat -ascii Inf_NOC_Na_fix
save output3\Inf_NOC_Nb.dat -ascii Inf_NOC_Nb_fix
save output3\Inf_NOC_Nele.dat -ascii Inf_NOC_Nele_fix
save output3\NOC_pm.dat -ascii NOC_pm_fix
save output3\XofN_pm.dat -ascii XofN_pm_fix
% dir = 'output4';
% mkdir(dir)
%save output4\Inf_NOC_Nele.txt -ascii Inf_NOC_Nele_fix
% writematrix(Inf_NOC_Nele2,'output4\element_set2.txt')%abaqusPML单元
ret2 = 'finish2';
