function [SE_pm] = PML_K(L,Lp,Fre,cr,G)
syms x
L=11; Lp=10;Fre=3;
f1x=10*(x-L)/Lp;
f2y=0;
ks=Fre/cr;
lamd1=(1+f1x/ks)-1i*f1x/ks;
lamd2=(1+f2y/ks)-1i*f2y/ks;
fm=lamd1*lamd2;
A_prime=[lamd2,0;0,lamd1];
Aprime=[1/lamd1,0;0,1/lamd2];
D_t=(A_prime*Aprime);
Gs=dlmread('Mate_S.txt',',',[1,0,1,0]);
const_matrices=Gs*D_t;

global XofNs NOC_Solid solid_layer1  Gs rs
Y1 = XofNs(NOC_Solid(n, 1), 2);
Y4 = XofNs(NOC_Solid(n, 4), 2);
% if y_a<=solid_layer1
%     lamds=dlmread('material_constant_Solid.txt',',',[1,0,1,0]);
%     Gs=dlmread('material_constant_Solid.txt',',',[1,1,1,1]);
% elseif y_a>solid_layer1
%     lamds=dlmread('material_constant_Solid.txt',',',[2,0,2,0]);
% end
SE = zeros (4,4);
XNI(1, 1) = -0.57735026919;
XNI(2, 1) = -0.57735026919;
XNI(1, 2) =  0.57735026919;
XNI(2, 2) = -0.57735026919;
XNI(1, 3) =  0.57735026919;
XNI(2, 3) =  0.57735026919;
XNI(1, 4) = -0.57735026919;
XNI(2, 4) =  0.57735026919;
for p = 1 :4;
    % nodal coordinates
    X1 = XofNs(NOC_Solid(n, 1), 1);
    Y1 = XofNs(NOC_Solid(n, 1), 2);
    X2 = XofNs(NOC_Solid(n, 2), 1);
    Y2 = XofNs(NOC_Solid(n, 2), 2);
    X3 = XofNs(NOC_Solid(n, 3), 1);
    Y3 = XofNs(NOC_Solid(n, 3), 2);
    X4 = XofNs(NOC_Solid(n, 4), 1);
    Y4 = XofNs(NOC_Solid(n, 4), 2);
    %%%%%%%%%%%integration points
    XI = XNI(1,p);
    ETA = XNI(2,p);
    TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4;
    TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4;
    TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4;
    TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4;
    DJ = TJ11 * TJ22 - TJ12 * TJ21;
    %  --- Inversion of JACOBIAN
    A(1, 1) = TJ22 / DJ;
    A(1, 2) = -TJ12 / DJ;
    A(2, 1) = -TJ21 / DJ;
    A(2, 2) = TJ11 / DJ;
    %  --- Local derivatives of shape function
    g(1, 1) = -(1 - ETA) / 4;
    g(2, 1) = -(1 - XI) / 4;
    g(1, 2) = (1 - ETA) / 4;
    g(2, 2) = -(1 + XI) / 4;
    g(1, 3) = (1 + ETA) / 4;
    g(2, 3) = (1 + XI) / 4;
    g(1, 4) = -(1 + ETA) / 4;
    g(2, 4) = (1 - XI) / 4;
    %  --- B = A * g
    B = zeros(2,4);
    for i = 1 : 2;
        for j = 1 : 4;
            c1 = 0;
            for k = 1 : 2;
                c1 = c1 + A(i, k) * g(k, j);
            end
            B(i , j) = c1;
        end
    end
    %  --- BB = pose the B in the strainBB
    % Element Stiffness Matrix  SE
    ShapeFunM(1,1)=0.25*(1-XI)*(1-ETA);
    ShapeFunM(1,2)=0.25*(1+XI)*(1-ETA);
    ShapeFunM(1,3)=0.25*(1+XI)*(1+ETA);
    ShapeFunM(1,4)=0.25*(1-XI)*(1+ETA);

    x=ShapeFunM*[X1;X2;X3;X4];
    C11=eval(const_matrices(x));
    SE=SE+B'*C11* B;
end
end


