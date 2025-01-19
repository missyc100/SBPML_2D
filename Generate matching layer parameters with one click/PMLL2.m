function [MU,CU,KU, CUS,KUS,MS,CS,KS,CSU,KSU] ...
    = PMLL2(n,XofN_pm,Inf_NOC_Na,Inf_NOC_Nb,ele_first,Infinite_Fault,R,m,ns,C_ref,b_ele,Lpml,L1,L2,PHY_LPML,AR1)

%%  没法进行积分
% m=2;
% Lpml=0.2;
% ns=1;
% C_ref=gama0*2*sqrt(Gs/rs)/(log(1/abs(R)));
% if ismember(n,Infinite_Fault)==1;%
%     lamds=dlmread('Mate_fault.txt',',',[1,0,1,0]);
%     Gs =dlmread('Mate_fault.txt',',',[1,1,1,1]);
%     rs =dlmread('Mate_fault.txt',',',[1,2,1,2]);
% else
%     lamds=dlmread('Mate_S.txt',',',[1,0,1,0]);
%     Gs =dlmread('Mate_S.txt',',',[1,1,1,1]);
%     rs =dlmread('Mate_S.txt',',',[1,2,1,2]);
% end
v=0.3;
Gs=1/2;
Et=2*Gs*(1+v);
lamds=2*Gs*v/(1-2*v);
rs=0.001;
% E=Gs.*(3*lamds+2*Gs)./(lamds+Gs);
%
% sqrt(Gs/rs)
% Gs =dlmread('Mate_S.txt',',',[1,0,1,0]);
% rs =dlmread('Mate_S.txt',',',[1,1,1,1]);
% vt=lamds./(2*(lamds+Gs));
% Et=Gs.*(3*lamds+2*Gs)./(lamds+Gs);
vt=v;
Dt=Et/((1+vt)*(1-2*vt))*[(1-vt),vt,0;vt,1-vt,0;0,0,0.5-vt];
Elastic=inv(Dt);


% Coo=[xbsa11,xbsa11,xbsa22,xbsa22;
%     -1,1 ,1,-1];
xbsa11=ele_first(n,1);
xbsa22=ele_first(n,2);
% Coo=[xbsa11,xbsa11,xbsa22,xbsa22;
%     -1,1 ,1,-1];
Coo=[xbsa11,xbsa11,xbsa22,xbsa22;
    1,-1 ,-1,1];


X1 = Coo(1,1);Y1 = Coo(2,1);
X2 = Coo(1,2);Y2 = Coo(2,2);
X3 = Coo(1,3);Y3 = Coo(2,3);
X4 = Coo(1,4);Y4 = Coo(2,4);



XNI(1, 1) = -0.57735026919;
XNI(2, 1) = -0.57735026919;
XNI(1, 2) =  0.57735026919;
XNI(2, 2) = -0.57735026919;
XNI(1, 3) =  0.57735026919;
XNI(2, 3) =  0.57735026919;
XNI(1, 4) = -0.57735026919;
XNI(2, 4) =  0.57735026919;

XB = XofN_pm(Inf_NOC_Nb(n, 1:2).', 1);
YB = XofN_pm(Inf_NOC_Nb(n, 1:2).', 2);
%XB=[0,1;1,0]*XB;
%YB=[0,1;1,0]*YB;
AR2 = [-XB(1)+XB(2);-YB(1)+YB(2)];
AR1(3) = [];
flag1 = dot(AR1,AR2);
if flag1<0
    XB=[0,1;1,0]*XB;
    YB=[0,1;1,0]*YB;
end

X0 = XofN_pm(Inf_NOC_Na(n, 1:2).', 1);
Y0 = XofN_pm(Inf_NOC_Na(n, 1:2).', 2);
AR3 = [-X0(1)+X0(2);-Y0(1)+Y0(2)];
flag2 = dot(AR1,AR3);
if flag2<0
    X0=[0,1;1,0]*X0;
    Y0=[0,1;1,0]*Y0;
end



% beta0=[120];
% alfa0=[0.75];
% beta0=20*100/3;% (m+1)*C_ref/(2*(Lpml*3))*log(1/abs(R));
% alfa0=1;%(m+1)*b_ele/(2*(Lpml*3))*log(1/abs(R));
beta0= (m+1)*C_ref/(2*(PHY_LPML))*log(1/abs(R));
alfa0=(m+1)*b_ele/(2*(PHY_LPML))*log(1/abs(R));


% Shape(1,1)=(xbsa2(1,1)-xbsa)*(1-eta)/(2*dxbas);
% Shape(2,1)=(xbsa2(1,1)-xbsa)*(1+eta)/(2*dxbas);
% Shape(3,1)=(xbsa-xbsa1(1,1))*(1-eta)/(2*dxbas);
% Shape(4,1)=(xbsa-xbsa1(1,1))*(1+eta)/(2*dxbas);

% B_xbsa_eta=[diff(Shape.',xbsa);diff(Shape.',eta)];


% double(subs(Yakobi_det,[eta,xbsa],[1,2]))
% double(Yakobi_det)
% x22=X0+xbsa22*((XB-X0));
% y22=Y0+xbsa22*((YB-Y0));
% x11=X0+xbsa11*((XB-X0));
% y11=Y0+xbsa11*((YB-Y0));

% coorxy=[x11.',x22.';y11.',y22.'];


% Inv_Yakobi=A_eta./Yakobi_det+xbsa*B_eta./Yakobi_det;
Mu=zeros(8); Cu=zeros(8); Ku=zeros(8);
CuS=zeros(8,12);
KuS=zeros(8,12);
Ks=zeros(12); Ms=zeros(12); Cs=zeros(12);


Zer2=zeros(2,1);

for p = 1 :4;
    % nodal coordinates
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
    Bbb = A * g;
    
    
    ShapeFunM(1,1)=0.25*(1-XI)*(1-ETA);
    ShapeFunM(1,2)=0.25*(1+XI)*(1-ETA);
    ShapeFunM(1,3)=0.25*(1+XI)*(1+ETA);
    ShapeFunM(1,4)=0.25*(1-XI)*(1+ETA);
    
    Bb1=[Bbb(:,1),Zer2,Bbb(:,2),Zer2,Bbb(:,3),Zer2,Bbb(:,4),Zer2];
    Bb2=[Zer2,Bbb(:,1),Zer2,Bbb(:,2),Zer2,Bbb(:,3),Zer2,Bbb(:,4)];
    
    NN1=[ShapeFunM(1,1),0,ShapeFunM(1,2),0,ShapeFunM(1,3),0,ShapeFunM(1,4),0];
    NN2=[0, ShapeFunM(1,1),0,ShapeFunM(1,2),0,ShapeFunM(1,3),0,ShapeFunM(1,4)];
    N12=[NN1;NN2];
    NN3=[ShapeFunM(1,1)*eye(3),ShapeFunM(1,2)*eye(3),ShapeFunM(1,3)*eye(3),ShapeFunM(1,4)*eye(3)];
    
    %     a111=(double(subs(ShapeFunM,[xbsa eta],[xbsa1,eta1])));
    
    
    
    
    
    xbsa1=ShapeFunM*Coo(1,:).';
    eta1=ShapeFunM*Coo(2,:).';
    eta=eta1;
    xbsa=xbsa1;
    
    N1=(1+eta)/2; N2=(1-eta)/2;
    Neta=[N1,N2];
    diff_eta=[1/2,-1/2];
    j11=(Neta*(XB-X0)*(diff_eta*Y0));
    j12=(Neta*(XB-X0)*(diff_eta*YB));
    j21=(Neta*(YB-Y0)*(diff_eta*X0));
    j22=(Neta*(YB-Y0)*(diff_eta*XB));
    
    a_eta=j11-j21;
    b_eta=-a_eta+j12-j22;
    % Yakobi_det=a_eta+b_eta*xbsa;
    
    A_eta=[diff_eta*Y0,-Neta*(YB-Y0);
        -diff_eta*X0,Neta*(XB-X0)];
    B_eta=[diff_eta*(YB-Y0),0;
        -diff_eta*(XB-X0),0];
    
    
    XX0=0;
    alfa_xbsa=1+alfa0*((xbsa-XX0)*ns/Lpml)^m;
    beta_xbsa=beta0*((xbsa-XX0)*ns/Lpml)^m;
    alfa_xb11=xbsa+alfa0*(ns/Lpml)^m*(xbsa-XX0)^(m+1)/(m+1);
    beta_xb11=beta0*(ns/Lpml)^m*(xbsa-XX0)^(m+1)/(m+1);
    
    % xbsa  xbsa
    % alfa_xbsa=1;
    % beta_xbsa=0;
    %
    % alfa_xb11=xbsa;
    % beta_xb11=0;
    
    Ea=[1,0;0,alfa_xbsa];
    Eb=[0,0;0,beta_xbsa];
    Fa=alfa_xb11*[1,0;0,1];
    Fb=beta_xb11*[1,0;0,1];
    
    a11=a_eta*alfa_xbsa+b_eta*alfa_xb11*alfa_xbsa;
    a22=a_eta*beta_xbsa+alfa_xb11*beta_xbsa*b_eta+b_eta*beta_xb11*alfa_xbsa;
    a33=b_eta*beta_xb11*beta_xbsa;
    
    
    a111=a11;%(double(subs(a11,[xbsa eta],[xbsa1,eta1])));
    a222=a22;%(double(subs(a22,[xbsa eta],[xbsa1,eta1])));
    a333=a33;%(double(subs(a33,[xbsa eta],[xbsa1,eta1])));
    %     B_xbsa_eta1=(double(subs(B_xbsa_eta,[xbsa eta],[xbsa1,eta1])));
    Ea1=Ea;%(double(subs(Ea,[xbsa eta],[xbsa1,eta1])));
    Eb1=Eb;%(double(subs(Eb,[xbsa eta],[xbsa1,eta1])));
    Fa1=Fa;%(double(subs(Fa,[xbsa eta],[xbsa1,eta1])));
    Fb1=Fb;%(double(subs(Fb,[xbsa eta],[xbsa1,eta1])));
    %     Yakob=(double(subs(Inv_Yakobi,[xbsa eta],[xbsa1,eta1])));
    
    %         Bbb1 = Yakob*A * g;
    %
    %     k=1:4;
    %     BBb(1,2*k-1)=Bbb1(1,k);
    %     BBb(3,2*k-1)=Bbb1(2,k);
    %     BBb(2,2*k)=Bbb1(2,k);
    %     BBb(3,2*k)=Bbb1(1,k);
    
    A_eta1=A_eta;%(double(subs(A_eta,[xbsa eta],[xbsa1,eta1])));
    B_eta1=B_eta;%(double(subs(B_eta,[xbsa eta],[xbsa1,eta1])));
    A1=A_eta1(1,:); A2=A_eta1(2,:); B1=B_eta1(1,:);B2=B_eta1(2,:);
    % B2*
    % Fb*B_eta*Eb
    CuS=CuS+[Bb1.'*[(A1*Ea1+B1*Fa1*Ea1).'],Bb2.'*[(A1*Ea1+B1*Fa1*Ea1).']]*L1.'*NN3*DJ+...
        [Bb1.'*[(A2*Ea1+B2*Fa1*Ea1).'],Bb2.'*[(A2*Ea1+B2*Fa1*Ea1).']]*L2.'*NN3*DJ;
    
    KuS=KuS+[Bb1.'*[(A1*Eb1+B1*Fa1*Eb1+B1*Fb1*Ea1).'],Bb2.'*[(A1*Eb1+B1*Fa1*Eb1+B1*Fb1*Ea1).']]*L1.'*NN3*DJ+...
        [Bb1.'*[(A2*Eb1+B2*Fa1*Eb1+B2*Fb1*Ea1).'],Bb2.'*[(A2*Eb1+B2*Fa1*Eb1+B2*Fb1*Ea1).']]*L2.'*NN3*DJ;
    
    
    Mu=Mu+rs*N12.'*a111*N12*DJ;
    Cu=Cu+rs*N12.'*a222*N12*DJ;
    Ku=Ku+rs*N12.'*a333*N12*DJ;
    
    % CuS1=CuS1+Bbb.'*Ea1.'*A_eta1.'*N2*DJ;
    % CuS2=CuS2+Bbb.'*Ea1.'*Fa1*B_eta1.'*N2*DJ;
    %
    %
    % KuS1=KuS1+Bbb.'*Eb1.'*A_eta1.'*N2*DJ;
    % KuS2=KuS2+Bbb.'*Ea1.'*Fb1*B_eta1.'*N2*DJ;
    % KuS3=KuS3+Bbb.'*Eb1.'*Fa1*B_eta1.'*N2*DJ;
    
    Ms=Ms+(NN3.')*a111*Elastic*NN3*DJ;
    Cs=Cs+(NN3.')*a222*Elastic*NN3*DJ;
    Ks=Ks+(NN3.')*a333*Elastic*NN3*DJ;
    
    %     Kkk=Kkk+BBb.'*Dt*BBb*DJ;

end

% Ms=1/Gs*double(int(int((N2.')*a11*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
% Cs=1/Gs*double(int(int((N2.')*a22*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
% Ks=1/Gs*double(int(int((N2.')*a33*N2, xbsa,xbsa11,xbsa22),eta,-1,1));


% CuS1=double(   int(  int( B_xbsa_eta.'*Ea.'*A_eta.'*N2,xbsa,xbsa11,xbsa22),eta,-1,1));
% CuS2=double(   int(  int( B_xbsa_eta.'*Ea.'*Fa*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
% CuS=CuS1+CuS2;
%
% % KuS1=double(   int(  int( B_xbsa_eta.'*Eb.'*A_eta.'*N2,xbsa,xbsa11,xbsa22),eta,-1,1));
% % KuS2=double(   int(  int( B_xbsa_eta.'*Ea.'*Fb*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
% % KuS3=double(   int(  int( B_xbsa_eta.'*Eb.'*Fa*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
% KuS=KuS1+KuS2+KuS3;

%
%
%
%
% % TuS=zeros(4,8);%double(   int(  int( B_xbsa_eta.'*Eb.'*Fb*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
CSu=CuS.';
KSu=KuS.';

Plusor=diag(Mu);
if all(Plusor(:)<0)==1
    Mu=-Mu;
    Cu=-Cu;
    Ku=-Ku;
    CuS=-CuS;
    KuS=-KuS;
    %  Kkk=-Kkk;
end

Plusor2=diag(Ms);
if all(Plusor2(:)<0)==1
    Ms=-Ms;
    Cs=-Cs;
    Ks=-Ks;
    CSu=-CSu;
    KSu=-KSu;
    %  Kkk=-Kkk;
end
% sum()

% CSu=CuS.';
% %
% % KSu2=double(   int(  int( N2.'*Fb*B_eta*Ea*B_xbsa_eta, xbsa,xbsa11,xbsa22),eta,-1,1));
%
% KSu=KuS.';


% DeDe=[2*nodXY11(:,1)-1;nodXY22(:,1)];
% de1=[nodXY11(:,1);nodXY222(:,1)];
% dd=[[1:4]',de1];
% de11=sortrows(dd,2);
% de=de11(:,1);
% De=zeros(8,1);
% De(1:2:end,1)=2*de-1;
% De(2:2:end,1)=2*de;
% De1=zeros(12,1);
% De1(1:3:end,1)=3*de-2;
% De1(2:3:end,1)=3*de-1;
% De1(3:3:end,1)=3*de;

MU=Mu;
CU=Cu;
KU=Ku;
CUS=CuS;
KUS=KuS;

MS=Ms;
CS=Cs;
KS=Ks;
CSU=CSu;
KSU=KSu;
% kkk=Kkk(De,De);
end

%     a_coefficient=(eval(subs(a1,[x y],[xx0,yy0])));
%     b_coefficient=(eval(subs(b1,[x y],[xx0,yy0])));
%     c_coefficient=(eval(subs(C1,[x y],[xx0,yy0])));

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
% Element Stiffness Matrix  SE
%     for i = 1 : 4;
%         for j = 1 : 4;
%             c1 = 0;
%             for k = 1 : 2;
%                 c1 = c1 + B(k, i)* B(k, j) * DJ;
%             end
%             SE(i, j) = c1+SE(i, j);
%         end
%     end
% end
% SE= Gs*SE;
% end