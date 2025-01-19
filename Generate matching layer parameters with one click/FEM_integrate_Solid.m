function [M_r,K_r] = FEM_integrate_Solid(NOC_INN,XofN_R,Finite_Fault,de2)
NE_INN=size(NOC_INN,1); NNs=size(XofN_R,1);
Ns=2*NNs;
AKE_s = zeros(Ns,Ns); %ASE表示 总刚矩阵
AMS_s= zeros(Ns,Ns); %AME表示 总刚矩阵
%%%%%%初始矩阵值
for n = 1 :NE_INN;
    SE=Ks(n,NOC_INN,XofN_R,Finite_Fault) ;
    ME=Ms(n,NOC_INN,XofN_R,Finite_Fault);
    %------形成总刚矩阵 K  M
    node = zeros(1,8);
    node(1:2:end) = 2*NOC_INN(n,:)-1;
    node(2:2:end) = 2*NOC_INN(n,:);

    AKE_s(node , node) = AKE_s(node , node) + SE;  %广播变量
    AMS_s(node , node) = AMS_s(node , node) + ME;
end



K_r=AKE_s(de2,de2);
M_r=AMS_s(de2,de2);


end
