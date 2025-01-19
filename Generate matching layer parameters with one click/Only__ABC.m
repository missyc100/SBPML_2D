% 两个Parts，包括 内域 与 PML 域
function ret1 = Only__ABC(NNr,NEr,NNpm,NEpm,rows)
% close all;clear all;
% clc
%% 内域
% NNr= 187; NEr=160;
NOC_R=dlmread('output1/element.txt',',',[0,1,NEr-1,4]);
XofN_R=dlmread('output1/node.txt',',',[0,1,NNr-1,2]);
% mesh=0.02m
N_F_Inter= importdata('output1/finite_inter.txt'); %% 有限域与无限域相交结点号
N_F_Inter= N_F_Inter(find( ~ isnan(N_F_Inter)));
%N_F_Inter=[1:201:2011].';
Infor_F_Inter=[N_F_Inter,XofN_R(N_F_Inter,:)];
Infor_F_I=sortrows(Infor_F_Inter, [2,3]);
NR=setdiff([1:NNr],Infor_F_I(:,1)).';
%% 外荷载
t1 = 1; c = 0.2;
a = @(t) ((t-t1)/c/t1).^2;
p = @(t) (1-a(t)).*exp(-0.5*a(t));
Tp=10;ntime=2001;
t=linspace(0,Tp,ntime);
loading=0.1*p(t);
dt=Tp/(ntime-1);
time=t.';
F0=loading;
%频谱
% F0 =fft(loading);  
% omiga=[0:2*pi/(Tp):2*pi/(dt)].';
% plot(omiga(1:end,:)./(2*pi),abs(fft(loading')));
Finite_Fault= importdata('Fi_Fault.txt'); % 此处为单元号
Finite_Fault= Finite_Fault(find( ~ isnan(Finite_Fault))); 
% 计算内域矩阵
de=[NR;Infor_F_I(:,1)];
de2=zeros(2*size(de,1),1);
de2(1:2:end-1,1)=2*de-1;
de2(2:2:end,1)=2*de;
% [Fi_MM,Fi_KK] = FEM_integrate_Solid(NOC_R,XofN_R,Finite_Fault,de2);  %J denotes the lumped
% ele_vec2=diff(a_B2_inf(:,2:2));
%% 加载
% For_De= importdata('output1/Force.txt');
% For_De= For_De(find( ~ isnan(For_De)));
% % For_De=[1:11:177].';
% [CoorXY,pressure] = enforce_lump(F0,XofN_R,For_De);
% CoorXY = CoorXY+1;%转为竖直方向
% [isornot,Force_N]=ismember(CoorXY,de2);

%% PML 层信息
% NNpm=259;
% NEpm=216 ;
NOC_pm=dlmread('output1/element_pml.txt',',',[0,1,NEpm-1,4]);
XofN_pm=dlmread('output1/node_pml.txt',',',[0,1,NNpm-1,2]);
N_I_Inter= importdata('output1/infinite_inter.txt'); %% 无限域与有限域相交结点号
N_I_Inter= N_I_Inter(find( ~ isnan(N_I_Inter)));
% N_I_Inter= [11:11:121].';  
Infor_I_Inter=[N_I_Inter,XofN_pm(N_I_Inter,:)];
Infor_I_I=sortrows(Infor_I_Inter, [2,3]);
N_PML=setdiff([1:NNpm],Infor_I_I(:,1)).';
De_PML=[Infor_I_I(:,1);N_PML];
Infinite_Fault= importdata('Infi_Fault.txt');
Infinite_Fault= Infinite_Fault(find( ~ isnan(Infinite_Fault)));
%                                               NOC_Solid,XofNs,Gs,rs,I_I,o_pml
% 边界层
Set_b= importdata('output1/Artificial_boundary.txt');
Set_b= Set_b(find( ~ isnan(Set_b)));
% Set_b=[1:11:111].'; % 固定结点号
Set_b=sortrows(Set_b); 
rows=rows+1;  % 层数+1；
[ Inf_NOC_Na, Inf_NOC_Nb,Inf_NOC_Nele,ele_first] = SBPML_NOC( NOC_pm,XofN_pm,rows,Set_b);
Nb1=Inf_NOC_Nb(:,1);Inf_NOC_Nb(:,1)=Inf_NOC_Nb(:,2);  
Inf_NOC_Nb(:,2)=Nb1;
Noc3=Inf_NOC_Nele(:,3);Inf_NOC_Nele(:,3)=Inf_NOC_Nele(:,4);
Inf_NOC_Nele(:,4)=Noc3;
% 刚度矩阵计算
[Inf_NOC_Nele2] = Assembly(Inf_NOC_Na,Inf_NOC_Nb,Inf_NOC_Nele,ele_first,XofN_pm,...
    Infinite_Fault,De_PML,Set_b);
% SN=2*size(NR,1)+5*(size(N_PML,1)+size(Infor_I_I,1));
%% 有限域无限域集成
% [MM,CC,KK]=Assembly_Finite_infinite(Fi_MM,Fi_KK,MM1,CC1,KK1,SN,NR);
% Real_Fix=Fix+2*size(NR,1);
%% Newmark 计算
%  [ UU,Dis ] = newmark( MM,CC,KK,ntime,dt,NNr,Real_Fix,Force_N,pressure);
%% 数据提取
% xxx=[For_De;Infor_F_I(:,1)];
% [isornot,Id]=ismember(xxx,[NR;Infor_F_I(:,1)]);%
% de=[2*Id];
% Uy_f=Dis(:,de);
% result=[t(1,:).',Uy_f(:,:)];
% figure
% plot(result(:,1),result(:,2))

%% 对应位置后处理
%writematrix(NOC_R,'A.txt')

dir = 'output2';
mkdir(dir)
save output2\ele_first.dat -ascii ele_first
save output2\Inf_NOC_Na.dat -ascii Inf_NOC_Na
save output2\Inf_NOC_Nb.dat -ascii Inf_NOC_Nb
save output2\Inf_NOC_Nele.dat -ascii Inf_NOC_Nele2
save output2\NOC_pm.dat -ascii NOC_pm
save output2\XofN_pm.dat -ascii XofN_pm

ret1 = 'finish1';



