function SEE=Ks1(n)
global XofNs NOC_Solid solid_layer1  Gs rs
Y1 = XofNs(NOC_Solid(n, 1), 2);
Y4 = XofNs(NOC_Solid(n, 4), 2);
Gs=dlmread('Mate_S.txt',',',[1,0,1,0]);
SEE = zeros (4,4);
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
    A=[TJ22 ,-TJ12;-TJ21,TJ11]/DJ;
    g=[ -(1 - ETA),(1 - ETA),(1 + ETA),-(1 + ETA);-(1 - XI),-(1 + XI), (1 + XI),(1 - XI)]./4;
    B= A * g;
   SEE=SEE+Gs*B'*B * DJ;
end
end



