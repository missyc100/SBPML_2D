function [SE_PML,ME_PML] = PML_M(n,L,Lp,Fre,cs,Gs,NOC_PML,rs)
global XofNs
syms x1 x2
f1x=20*(x1-L)/Lp;
Q=[sqrt(2)/2,-sqrt(2)/2;sqrt(2)/2,sqrt(2)/2];
inv_Q=inv(Q);
f2y=0;
ks=Fre/cs;
lamd1=(1+f1x/ks)-1i*f1x/ks;
lamd2=(1+f2y/ks)-1i*f2y/ks;
fm=lamd1*lamd2;
A_prime=[lamd2,0;0,lamd1];
Aprime=[1/lamd1,0;0,1/lamd2];
D_t1=Q*A_prime*Aprime*Q';
const_matrices=Gs*D_t1;
SE_PML = zeros (4,4);
ME_PML = zeros (4,4);

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
    X1 = XofNs(NOC_PML(n, 1), 1);
    Y1 = XofNs(NOC_PML(n, 1), 2);
    X2 = XofNs(NOC_PML(n, 2), 1);
    Y2 = XofNs(NOC_PML(n, 2), 2);
    X3 = XofNs(NOC_PML(n, 3), 1);
    Y3 = XofNs(NOC_PML(n, 3), 2);
    X4 = XofNs(NOC_PML(n, 4), 1);
    Y4 = XofNs(NOC_PML(n, 4), 2);
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
    y=ShapeFunM*[Y1;Y2;Y3;Y4];
    prime_x=inv_Q*[x;y];
    x1=prime_x(1,1);
    C11=eval(real(const_matrices));
    C11_I=eval(imag(const_matrices));
    CCC=C11+1i*C11_I;
    Fmm=eval(real(fm));
    Fmm_I=eval(imag(fm));
    Fm1=Fmm+1i*Fmm_I;
    SE_PML=SE_PML+B'*CCC* B*DJ;
    ME_PML=ME_PML+rs*Fm1*ShapeFunM'*ShapeFunM*DJ;
end
end


