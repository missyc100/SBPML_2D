function [Arr,Ars,Ass,MM] = Ele_matrix(n,B_infLR,nxLR)
Y=B_infLR(:,3);
Y1= Y(n);Y2= Y(n+1);
nx=nxLR(n,1);ny=nxLR(n,2);
Lr=-1/ny*[-ny,nx]';
Ls=1/ny*[0,1]';
G=dlmread('Mate_S.txt',',',[n,0,n,0]);
den=dlmread('Mate_S.txt',',',[n,1,n,1]);
syms s
s1=Y1/ny; s2=Y2/ny;
N1=(s2-s)/(s2-s1); N2=(s-s1)/(s2-s1);
Shape=[N1,N2];
D_Shape_s=diff(Shape,s);
Arr=eval(ny*int(Shape'*(G*Lr'*Lr)*Shape, s1,s2));
Ars=eval(ny*int(Shape'*(G*Lr'*Ls)*D_Shape_s, s1,s2));
Ass=eval(ny*int(D_Shape_s'*(G*Ls'*Ls)*D_Shape_s, s1,s2));
MM=eval(ny*den*int(Shape'*Shape, s1,s2));
end

