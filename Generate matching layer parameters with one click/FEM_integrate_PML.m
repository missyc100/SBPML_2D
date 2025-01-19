function [UU] = FEM_integrate_PML(F0,NOC_PML,Fre,cs,Gs,rs,Fi_MM,Fi_KK)
global   NNs XofNs
NE_PML=size(NOC_PML,1);
AKE_pml = zeros(NNs,NNs); % ASE��ʾ �ܸվ���
AMS_pml= zeros(NNs,NNs);  % AME��ʾ �ܸվ���
L=1.5;Lp=0.5;
%%%%%%��ʼ����ֵ
for n = 1 :NE_PML;
    [SE_pml,ME_pml]=PML_M(n,L,Lp,Fre,cs,Gs,NOC_PML,rs);
    %------�γ��ܸվ��� K  M
    node = zeros(1,4);
    node = NOC_PML(n,:);
    AKE_pml(node , node) = AKE_pml(node , node) + SE_pml;  %�㲥����
    AMS_pml(node , node) = AMS_pml(node , node) + ME_pml;
end
% �ܺɵ�
[ro1,co]=find(XofNs(:,2)==0 & XofNs(:,1)==0);
% �̶��߽粿��
[ro,co]=find(XofNs(:,2)==max(XofNs(:,2)));

[fix_ro,co]=find(XofNs(:,1)==max(XofNs(:,1)));
ro=unique([ro;fix_ro]);
Fc=zeros(size(XofNs,1),1);


Fc(ro1,1)=F0;
Keq=-Fre^2*(AMS_pml+Fi_MM)+AKE_pml+Fi_KK;
Keq(ro,:)=0;
Keq(:,ro)=0;

for jj=1:length(ro)
    Keq(ro(jj),ro(jj))=1;
end
Fc(ro,1)=0;
UU=sparse(Keq)\Fc;
end
