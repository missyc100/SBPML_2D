function [ UU ] = Solve_E( AMS_pml,AKE_pml,Fi_MM,Fi_KK,Fre,F0,XofNs)
[ro,co]=find(XofNs(:,2)==max(XofNs(:,2)));
[ro1,co]=find(XofNs(:,2)==0 & XofNs(:,1)==0);

Fc=zeros(size(XofNs,1),1);
Fc(ro1,1)=F0;
Keq=-Fre^2*(AMS_pml+Fi_MM)+AKE_pml+Fi_KK;
Keq(ro,ro)=eye(size(co,1));
UU=sparse(Keq)\Fc;

end

