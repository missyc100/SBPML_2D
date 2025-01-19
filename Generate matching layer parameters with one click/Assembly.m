function [NOC_Solid] = Assembly(Inf_NOC_Na,Inf_NOC_Nb,Inf_NOC_Nele,ele_first,XofN_pm,...
    Infinite_Fault,N_PML,Set_b)

% NOC_Solid=NOC_pm;
% XofNs=XofN_pm;
% Fixed= importdata('Fixed.txt');
% N_B= Fixed(find( ~ isnan(Fixed)));
% 
% N_boundary=[N_B,XofNs(N_B,:)];
% N_boundary=sortrows(N_boundary, [2,3]);





NE=size(Inf_NOC_Nele,1);
NS=2*size(XofN_pm,1);
NS3=3*size(XofN_pm,1);


%AMU =zeros(NS,NS);ACU =AMU; AKU =AMU;
%ACUS =zeros(NS,NS3);AKUS =ACUS; 
%AMS =zeros(NS3,NS3);ACS =AMS; AKS =AMS; 
%ACSU =ACUS.';AKSU =AKUS.'; 

% MM=zeros(3*NS,3*NS);
% KK=zeros(3*NS,3*NS);
% CC=zeros(3*NS,3*NS);

%Zeross=zeros(NS,NS3);
% AA_alfay =AMA;AA_alfax =AMA; AA_betay =AMA; AA_betax =AMA;
% AA_alfay1 =AMA;AA_alfax1 =AMA; AA_betay1 =AMA; AA_betax1 =AMA;
% Au_zeros = AMA;
% C_ref=gama0*2*sqrt(Gs/rs)/(log(1/abs(R)));
% node = zeros(1,4); 
node2 = zeros(1,8); node3 = zeros(1,12);

% 匹配层参数设置位置
R=1.0e-10;
gama0=4;
b_ele=0.33*10;  % ele=0.02
m=2;
Lpml=1;
PHY_LPML=2;
ns=1;
C_ref= 22.36;

L1=[1,0;0,0;0,1];
L2=[0,0;0,1;1,0];
NOC_Solid=Inf_NOC_Nele;
%%%%%%初始矩阵值
% mm=zeros(NE,3)
for n =1 :NE;
    %     tic
    n
    AAX1 = XofN_pm(NOC_Solid(n, 1), 1);
    AAY1 = XofN_pm(NOC_Solid(n, 1), 2);
    ABX1 = XofN_pm(NOC_Solid(n, 2), 1);
    ABY1 = XofN_pm(NOC_Solid(n, 2), 2);
    ACX1 = XofN_pm(NOC_Solid(n, 3), 1);
    ACY1 = XofN_pm(NOC_Solid(n, 3), 2);
    AR1 = [-AAX1+ABX1;-AAY1+ABY1;0];
    AR2 = [-ABX1+ACX1;-ABY1+ACY1;0];
    flag1 = cross(AR1,AR2);
    if flag1(3,1)<0
        NOC_Solid(n,:) = NOC_Solid(n,:)*[0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0];
        AR1=-AR1;
    end
%     [ MU,CU,KU, CUS,KUS,MS,CS,KS,CSU,KSU]=PMLL2(n,C_ref,b_ele,XofNs,NOC_Solid,Infinite_Fault,R,m,Lpml,ns,X_Spline,N_boundary,L1,L2,PHY_LPML);
%     toc
%     [ MU,CU,KU, CUS,KUS,MS,CS,KS,CSU,KSU]=PMLL2(n,XofN_pm,Inf_NOC_Na,Inf_NOC_Nb,ele_first,Infinite_Fault,R,m,ns,C_ref,b_ele,Lpml,L1,L2,PHY_LPML,AR1);

 
%     tic
%     [ MU1,CU1,KU1, CUS1,KUS1,MS1,CS1,KS1,CSU1,KSU1]=PMLL(n,C_ref,b_ele,XofNs,NOC_Solid,R,m,Lpml,ns,X_Spline,Infor_I_I);
%     toc
    %     n
%     mm(n,1)=dxbas; mm(n,2)=xbsa22; mm(n,3)=xbsa11;
% end
    %------形成总刚矩阵 K  M
    
%     node = NOC_Solid(n,:);
%     node2(1,1:2:end)=2*node-1;
%     node2(1,2:2:end)=2*node;
%     node3(1,1:3:end)=3*node-2;
%     node3(1,2:3:end)=3*node-1;
%     node3(1,3:3:end)=3*node;

    
    
%     AMU(node2 , node2) = AMU(node2 , node2) + MU;  %广播变量
%     ACU(node2 , node2) = ACU(node2 , node2) + CU;
%     AKU(node2 , node2) = AKU(node2 , node2) + KU;
%     
%     ACUS(node2 , node3) = ACUS(node2 , node3) + CUS;  %广播变量
%     AKUS(node2 , node3) = AKUS(node2 , node3) + KUS;
%     
%     ACSU(node3 , node2) = ACSU(node3 , node2) + CSU;  %广播变量
%     AKSU(node3 , node2) = AKSU(node3 , node2) + KSU;
%     
%     
%     AMS(node3 , node3) = AMS(node3 , node3) + MS;
%     ACS(node3 , node3) = ACS(node3 , node3) + CS;
%     AKS(node3 , node3) = AKS(node3 , node3) + KS;
%     AA_betax(node , node) = AA_betax(node , node) + A_betax;
%     
%     AA_alfay1(node , node) = AA_alfay1(node , node) + A_alfay1;
%     AA_alfax1(node , node) = AA_alfax1(node , node) + A_alfax1;
%     AA_betay1(node , node) = AA_betay1(node , node) + A_betay1;
%     AA_betax1(node , node) = AA_betax1(node , node) + A_betax1;
end
% 受荷点
% MM=[AMU,Zeross;
%     Zeross.',-AMS];
% clear AMU AMS
% CC=[ACU,ACUS,;
%     ACSU,-ACS];
% clear ACU ACUS ACSU ACS
% KK=[AKU,AKUS;
%     AKSU,-AKS];
% clear AKU AKUS AKSU AKS
% 
% 
% de1=[N_PML];
% de2=zeros(NS,1);
% de2(1:2:end-1,1)=2*de1-1;
% de2(2:2:end,1)=2*de1;
% 
% De=[de2;[NS+1:5*NS/2].'];
% Mpml=MM(De,De);
% Cpml=CC(De,De);
% Kpml=KK(De,De);
% 
% 
% % esp=1.0e-9;
% % X_max=max(XofNs(:,1));X_min=min(XofNs(:,1));
% % Y_max=max(XofNs(:,2));
% % P=[X_min,X_max,Y_max].';
% dbx=2*Set_b-1; dby=2*Set_b;
% db=[dbx;dby];
% [isornot,Fix11]=ismember(db,De);%

end
