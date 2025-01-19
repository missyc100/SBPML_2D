function [ Absorb_matrices ] =bottom_absorb(XofBottom,Ah)

XB=XofBottom(:,2);
NOC_B=zeros(size(XB,1)-1,2);
NOC_B(:,1)=[1:size(XB,1)-1]'; NOC_B(:,2)=[2:size(XB,1)]';
Absorb_matrices=zeros(size(XB,1));
for n = 1:size(NOC_B,1)
    s1=XB(n); s2=XB(n+1);
        syms s
    N1=(s2-s)/(s2-s1); N2=(s-s1)/(s2-s1);
    Shape=[N1,N2];
    BE=eval(int(Shape'*Ah*Shape, s1,s2));
    degree = zeros(1,2);
    degree=NOC_B(n,:);
    Absorb_matrices(degree,degree) =  Absorb_matrices(degree,degree) + BE;
end
end


