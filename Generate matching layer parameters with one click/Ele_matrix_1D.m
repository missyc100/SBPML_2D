function [A11,A22,MM] = Ele_matrix_1D(n,B_infLR)
Y=B_infLR(:,1);
y1= Y(n);y2= Y(n+1);
G=dlmread('Mate_S.txt',',',[n,0,n,0]);
den=dlmread('Mate_S.txt',',',[n,1,n,1]);
syms y
N1=(y2-y)/(y2-y1); N2=(y-y1)/(y2-y1);
Shape=[N1,N2];
D_Shape_s=diff(Shape,y);
A11=eval(int(Shape'*(G)*Shape, y1,y2));
A22=eval(int(D_Shape_s'*(G)*D_Shape_s, y1,y2));
MM=eval(den*int(Shape'*Shape, y1,y2));
end

