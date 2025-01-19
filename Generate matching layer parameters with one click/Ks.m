function SE=Ks(n,NOC_INN,XofN_R,Finite_Fault)
% NOC_R,XofN_R,
XofNA=XofN_R ;
NOC_Solid=NOC_INN;
% Y1 = XofNA(NOC_Solid(n, 1), 2);
% Y4 = XofNA(NOC_Solid(n, 4), 2);
% y_a=(Y1+Y4)/2;
% % if y_a<=solid_layer1
% %     lamds=dlmread('material_constant_Solid.txt',',',[1,0,1,0]);
% %     Gs=dlmread('material_constant_Solid.txt',',',[1,1,1,1]);
% % elseif y_a>solid_layer1
% %     lamds=dlmread('material_constant_Solid.txt',',',[2,0,2,0]);
% %     Gs=dlmread('material_constant_Solid.txt',',',[2,1,2,1]);
% % end
% if ismember(n,Finite_Fault)==1;%
%     lamds=dlmread('Mate_fault.txt',',',[1,0,1,0]);
%     Gs=dlmread('Mate_fault.txt',',',[1,1,1,1]);
% else
%     lamds=dlmread('Mate_S.txt',',',[1,0,1,0]);
%     Gs=dlmread('Mate_S.txt',',',[1,1,1,1]);
% end




v=0.3;
Gs=1/2;
E=2*Gs*(1+v);
lamds=2*Gs*v/(1-2*v);
% E=Gs.*(3*lamds+2*Gs)./(lamds+Gs);
D=E/((1+v)*(1-2*v))*[(1-v),v,0;v,1-v,0;0,0,0.5-v];
% Cp=sqrt((lamd+2*G)/1);
% Cs=sqrt((G)/1);



SE = zeros (8,8);
BB=zeros(3,8);
% BD=zeros(8,3);
% v=lamds./(2*(lamds+Gs));

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
    B = A * g;
%     
%     
%     B = zeros(2,4);
%     for i = 1 : 2;
%         for j = 1 : 4;
%             c1 = 0;
%             for k = 1 : 2;
%                 c1 = c1 + A(i, k) * g(k, j);
%             end
%             B(i , j) = c1;
%         end
%     end
    %  --- BB = pose the B in the strainBB
    k=1:4;
    BB(1,2*k-1)=B(1,k);
    BB(3,2*k-1)=B(2,k);
    BB(2,2*k)=B(2,k);
    BB(3,2*k)=B(1,k);
    
%     %  --- strain_BB'*D= B'D
%     for i=1:8
%         for j=1:3
%             c2=0;
%             for k=1:3
%                 c2=c2+BB(k,i)*D(k,j);
%             end
%             BD(i,j)=c2;
%         end
%     end
    BD=BB.'*D;
    %  --- SE=B'DB:the stiffness SE
%     for i=1:8
%         for j=1:8
%             c3=0;
%             for k=1:3
%                 c3=c3+BD(i,k)*BB(k,j)*DJ;
%             end
%             SE(i,j)=c3+SE(i,j);
%         end
%     end
    SE= SE+BD*BB*DJ;
%     SE2-SE
end
end