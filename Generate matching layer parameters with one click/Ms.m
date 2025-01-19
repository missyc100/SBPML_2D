function ME=Ms(n,NOC_INN,XofN_R,Finite_Fault)
XofNA=XofN_R;
NOC_Solid=NOC_INN;

% if ismember(n,Finite_Fault)==1;%
%     rs=dlmread('Mate_fault.txt',',',[1,2,1,2]);
% else
%     rs=dlmread('Mate_S.txt',',',[1,2,1,2]);
% end

% Y1 = XofNA(NOC_Solid(n, 1), 2);
% Y4 = XofNA(NOC_Solid(n, 4), 2);
% % y_a=(Y1+Y4)/2;
rs=0.001;
den2=rs;
% if y_a<=solid_layer1
%     den2 =dlmread('material_constant_Solid.txt',',',[1,2,1,2]);
% elseif y_a>solid_layer1
%     den2 =dlmread('material_constant_Solid.txt',',',[2,2,2,2]);
% end
ME = zeros (8,8);ME2 = zeros (8,8);

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
    X1 = XofNA(NOC_Solid(n, 1), 1);
    Y1 = XofNA(NOC_Solid(n, 1), 2);
    X2 = XofNA(NOC_Solid(n, 2), 1);
    Y2 = XofNA(NOC_Solid(n, 2), 2);
    X3 = XofNA(NOC_Solid(n, 3), 1);
    Y3 = XofNA(NOC_Solid(n, 3), 2);
    X4 = XofNA(NOC_Solid(n, 4), 1);
    Y4 = XofNA(NOC_Solid(n, 4), 2);
    %integration points
    XI = XNI(1,p);
    ETA = XNI(2,p);
    ShapeFunM(1,1)=0.25*(1-XI)*(1-ETA);
    ShapeFunM(1,2)=0.25*(1+XI)*(1-ETA);
    ShapeFunM(1,3)=0.25*(1+XI)*(1+ETA);
    ShapeFunM(1,4)=0.25*(1-XI)*(1+ETA);
    NN1=[ShapeFunM(1,1),0,ShapeFunM(1,2),0,ShapeFunM(1,3),0,ShapeFunM(1,4),0;
    0, ShapeFunM(1,1),0,ShapeFunM(1,2),0,ShapeFunM(1,3),0,ShapeFunM(1,4)];
%     N12=[NN1;NN2];

    TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4;
    TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4;
    TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4;
    TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4;
    DJ = TJ11 * TJ22 - TJ12 * TJ21;
%     for i=1:4
%         for j=1:4
%             ME(2*i-1,2*j-1)=ME(2*i-1,2*j-1)+den2(1,1)*ShapeFunM(1,i)*ShapeFunM(1,j)*DJ;
%             ME(2*i,2*j)=ME(2*i,2*j)+den2(1,1)*ShapeFunM(1,i)*ShapeFunM(1,j)*DJ;
%         end
%     end
    ME=ME+den2*NN1.'*NN1*DJ;

    
end
%     ME2-ME

end

     
           

