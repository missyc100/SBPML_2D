function [MU,CU,KU, CUS,KUS,MS,CS,KS,CSU,KSU] ...
    = PMLL(n,C_ref,b_ele,XofNs,NOC_Solid,R,m,Lpml,ns,X_Spline,Infor_I_I)
NOC_PML=NOC_Solid;

% m=2;
% Lpml=0.2;
% ns=1;
% C_ref=gama0*2*sqrt(Gs/rs)/(log(1/abs(R)));

Gs =dlmread('Mate_S.txt',',',[1,0,1,0]);
rs =dlmread('Mate_S.txt',',',[1,1,1,1]);


X1 = XofNs(NOC_PML(n, 1), 1);
Y1 = XofNs(NOC_PML(n, 1), 2);
X2 = XofNs(NOC_PML(n, 2), 1);
Y2 = XofNs(NOC_PML(n, 2), 2);
X3 = XofNs(NOC_PML(n, 3), 1);
Y3 = XofNs(NOC_PML(n, 3), 2);
X4 = XofNs(NOC_PML(n, 4), 1);
Y4 = XofNs(NOC_PML(n, 4), 2);
CoorXY=[X1,X2,X3,X4;
    Y1,Y2,Y3,Y4];
% ele_vec2=diff(a_B2_inf(:,2:3));
% unit_vec2=zeros(length(ele_vec2),2);
% for jj=1:length(ele_vec2)
%     unit_vec2(jj,:)=ele_vec2(jj,:)./ norm(ele_vec2(jj,:));
% end
%% 为寻找比例线坐标及其边界坐标定义单元直线
syms X Y
dx=diff([CoorXY(1,:),CoorXY(1,1)]);
dy=diff([CoorXY(2,:),CoorXY(2,1)]);
% 找出单元四条直线
Kk=dy./dx;
Kk=roundn(Kk,-3);
Kk(Kk(1,:) == -inf)=inf;
In_coor1=[1:4;CoorXY;Kk];
% In_coor1=roundn(In_coor,-5);
% find(In_coor1(4,:))
% A = [3 3 2 3]; %没有相同元素时，b为空
Kk_num = tabulate(In_coor1(4,:));
K_same = Kk_num(Kk_num(:,2) ~= 1 & Kk_num(:,2) ~= 0);
midy=0.6;

if size(K_same,1)>1
    K_same= sort(K_same);
    if mean(CoorXY(2,:))<midy  
        K_same=K_same(K_same~=0);
    else
        K_same=K_same(K_same==0);
    end
end

    
    
    
[a,co]=find(K_same==In_coor1(4,:));
Infor_coor=[In_coor1,In_coor1(:,1)];
% nodc=zeros(2,3);

num1=co(1):co(1)+1;
num2=co(2):co(2)+1;

    
X1=Infor_coor(2,num1).';
Y1=Infor_coor(3,num1).';
nodXY1=[num1',X1,Y1];
nodXY11=sortrows(nodXY1,[2,3]);

X2=Infor_coor(2,num2).';
Y2=Infor_coor(3,num2).';
nodXY2=[num2',X2,Y2];

nodXY22=sortrows(nodXY2,[2,3]);
nodXY22(nodXY22(:,1)>4)=1;

X1=nodXY11(:,2);
Y1=nodXY11(:,3);

X2=nodXY22(:,2);
Y2=nodXY22(:,3);





% num=unique(In_coor1(4,:));

% In_coor1=(sortrows(In_coor.',4)).';
% difKk1=diff(In_coor1(4,:));
% difKk1(abs(difKk1)<esp)=0;
% MindKk = find(difKk1==0);

% b = union(In_coor(4,:),[])
% In_coor=[1:4;CoorXY;Kk];

% index=find(abs(Kk(1,:))>1.0e10);
% Kk(1,index)=1.0e10;
% Kk1=setdiff(Kk,K_same)
[a,b]=find(Kk==K_same);
ro1=setdiff([1:4],b).';
F=Kk.*(X-CoorXY(1,:))+CoorXY(2,:)-Y;
if ismember(inf,Kk)==1;
    [ro,co]=find(Kk==inf);
    F(1,co)=(X-CoorXY(1,co));
end



esp=1.0e-5;
% signx=sign(mean(CoorXY(1,:)));
% midY=mean(CoorXY(2,:));
% % find(X_Spline(:,2)<)
Node_Spline=size(X_Spline,1);
XX=[];
YY=[];
XB1=[];
YB1=[];


% ro1
for ii=1:size(ro1,1)
    cc=ro1(ii,1);
    gg=0; mm=0;
    for jj=1:Node_Spline
        Value_spline=(double(subs(F(1,cc),[X, Y],[X_Spline(jj,2) ,X_Spline(jj,3) ])));
        Value_Bound=(double(subs(F(1,cc),[X, Y],[Infor_I_I(jj,2) ,Infor_I_I(jj,3) ])));
%         gg=0;
        if  abs(Value_spline)<1.0e-4
            gg=1+gg;
            XX(ii,gg)=X_Spline(jj,2);
            YY(ii,gg)=X_Spline(jj,3);
        end
%         gg=0;
        if  abs(Value_Bound)<1.0e-4
            mm=1+mm;
            XB1(ii,mm)=Infor_I_I(jj,2);
            YB1(ii,mm)=Infor_I_I(jj,3);
        end
    end
end

if size(XX,2)>1
    if mean(CoorXY(1,:))<0
        [ro,co]=(find(XX(:,:)<0));
        XX=XX(ro,co(1));
        YY=YY(ro,co(1));
        
        [ro,co]=(find(XB1(:,:)<0));
        XB1=XB1(ro,co(1));
        YB1=YB1(ro,co(1));
    else
        [ro,co]=(find(XX(:,:)>0));
        XX=abs(XX(ro,co(1)));
        YY=abs(YY(ro,co(1)));
        
        [ro,co]=(find(XB1(:,:)>0));
        XB1=abs (XB1(ro,co(1)));
        YB1=abs(YB1(ro,co(1)));
    end
end

XY0=[XX(:,:),YY(:,:)];
XY0=sortrows(XY0,[1,2]);
% YY=
XYB=[XB1(:,:),YB1(:,:)];
XYB=sortrows(XYB,[1,2]);



% if signx<0
%     [kx,cx] = find(XX<0);
%     [kbx,cbx] = find(XB1<0);
% else
%     [kx,cx] = find(XX>0);
%     [kbx,cbx] = find(XB1>0);
% end

% x01=XX(1,1); x02=XX(2,1);
% y01=YY(1,1); y02=YY(2,1);
X0=XY0(:,1);
Y0=XY0(:,2);

% xb1=XB1(kbx(1,1),cbx(1,1)); xb2=XB1(kbx(2,1),cbx(2,1));
% yb1=YB1(kbx(1,1),cbx(1,1)); yb2=YB1(kbx(2,1),cbx(2,1));
% XB=[xb1;xb2];
% YB=[yb1;yb2];
XB=XYB(:,1);
YB=XYB(:,2);
if abs (norm(XB-X0))<esp
%     X1=Y1; X2=Y2;
%     X0=Y0;
%     XB=YB;
    xbsa1=roundn((Y1-Y0)./(YB-Y0),-4);
    xbsa2=roundn((Y2-Y0)./(YB-Y0),-4);
else
    xbsa1=roundn((X1-X0)./(XB-X0),-4);
    xbsa2=roundn((X2-X0)./(XB-X0),-4);
end
while size(xbsa1,1)~=2 & size(xbsa2,1)~=2 
    break
end
if ismember (1,isnan(xbsa1))
    [a,b,c]=find(isnan(xbsa1));
    xbsa1( isnan(xbsa1))=[];
    xbsa1(a,b)=xbsa1;
    xbsa2( isnan(xbsa2))=[];
    xbsa2(a,b)=xbsa2;
    %     xbsa1()
end


if xbsa1(1,1)>xbsa2(1,1)
    aa=xbsa1;
    xbsa1=xbsa2;
    xbsa2= aa;
    
    bb=nodXY11;
    nodXY11=nodXY22;
    nodXY22= bb;
end
xbsa22=xbsa2(1,1);
xbsa11=xbsa1(1,1);

dxbas=xbsa2(1,1)-xbsa1(1,1);
while dxbas<0 
    break
end
while xbsa22>3 || xbsa11>3
    break
end

X11=nodXY11(:,2);
Y11=nodXY11(:,3);

X22=nodXY22(:,2);
Y22=nodXY22(:,3);



% In_coor1=sortrows(In_coor.',2).';


syms eta xbsa
N1=(1-eta)/2; N2=(1+eta)/2;
Neta=[N1,N2]; 
diff_eta=diff(Neta,eta);
j11=(Neta*(XB-X0)*(diff_eta*Y0));
j12=(Neta*(XB-X0)*(diff_eta*YB));
j21=(Neta*(YB-Y0)*(diff_eta*X0));
j22=(Neta*(YB-Y0)*(diff_eta*XB));

a_eta=j11-j21;
b_eta=-a_eta+j12-j22;
Yakobi_det=a_eta+b_eta*xbsa;
% double(subs(Yakobi_det,[eta,xbsa],[1,2]))
% double(Yakobi_det)
% x22=X0+xbsa22*((XB-X0));
% y22=Y0+xbsa22*((YB-Y0));
% x11=X0+xbsa11*((XB-X0));
% y11=Y0+xbsa11*((YB-Y0));

% coorxy=[x11.',x22.';y11.',y22.'];

A_eta=[diff_eta*Y0,-Neta*(YB-Y0);
    -diff_eta*X0,Neta*(XB-X0)];
B_eta=[diff_eta*(YB-Y0),0;
    -diff_eta*(XB-X0),0];
% Inv_Yakobi=A_eta./Yakobi_det+xbsa*B_eta./Yakobi_det;


ShapeFunM(1,1)=(xbsa2(1,1)-xbsa)*(1-eta)/(2*dxbas);
ShapeFunM(2,1)=(xbsa2(1,1)-xbsa)*(1+eta)/(2*dxbas);
ShapeFunM(3,1)=(xbsa-xbsa1(1,1))*(1-eta)/(2*dxbas);
ShapeFunM(4,1)=(xbsa-xbsa1(1,1))*(1+eta)/(2*dxbas);

N2=[ShapeFunM(1,1),0,ShapeFunM(2,1),0,ShapeFunM(3,1),0,ShapeFunM(4,1),0;
     0, ShapeFunM(1,1),0,ShapeFunM(2,1),0,ShapeFunM(3,1),0,ShapeFunM(4,1)];



beta0=(m+1)*C_ref/(2*Lpml)*log(1/abs(R));
alfa0=(m+1)*b_ele/(2*Lpml)*log(1/abs(R));

alfa_xbsa=1+alfa0*((xbsa-1)*ns/Lpml)^m;
beta_xbsa=beta0*((xbsa-1)*ns/Lpml)^m;
alfa_xb11=xbsa+alfa0*(ns/Lpml)^m*(xbsa-1)^(m+1)/(m+1);
beta_xb11=beta0*(ns/Lpml)^m*(xbsa-1)^(m+1)/(m+1);

% alfa_xbsa=1;
% beta_xbsa=0;
% 
% alfa_xb11=xbsa;
% beta_xb11=0;

Ea=[1,0;0,alfa_xbsa];
Eb=[0,0;0,beta_xbsa];
Fa=alfa_xb11*[1,0;0,1];
Fb=beta_xb11*[1,0;0,1];

% Eb.'*Fb*B_eta.'
% Fb*B_eta*Eb

a11=a_eta*alfa_xbsa+b_eta*alfa_xb11*alfa_xbsa;
a22=a_eta*beta_xbsa+alfa_xb11*beta_xbsa*b_eta+b_eta*beta_xb11*alfa_xbsa;
a33=b_eta*beta_xb11*beta_xbsa;

B_xbsa_eta=[diff(ShapeFunM.',xbsa);diff(ShapeFunM.',eta)];





Mu=rs*double(int(int((ShapeFunM)*a11*ShapeFunM.', xbsa,xbsa11,xbsa22),eta,-1,1));
Cu=rs*double(int(int((ShapeFunM)*a22*ShapeFunM.', xbsa,xbsa11,xbsa22),eta,-1,1));
Ku=rs*double(int(int((ShapeFunM)*a33*ShapeFunM.', xbsa,xbsa11,xbsa22),eta,-1,1));

CuS1=double(   int(  int( B_xbsa_eta.'*Ea.'*A_eta.'*N2,xbsa,xbsa11,xbsa22),eta,-1,1));
CuS2=double(   int(  int( B_xbsa_eta.'*Ea.'*Fa*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
CuS=CuS1+CuS2;

KuS1=double(   int(  int( B_xbsa_eta.'*Eb.'*A_eta.'*N2,xbsa,xbsa11,xbsa22),eta,-1,1));
KuS2=double(   int(  int( B_xbsa_eta.'*Ea.'*Fb*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
KuS3=double(   int(  int( B_xbsa_eta.'*Eb.'*Fa*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
KuS=KuS1+KuS2+KuS3;

CSu=CuS.';
KSu=KuS.';




% TuS=zeros(4,8);%double(   int(  int( B_xbsa_eta.'*Eb.'*Fb*B_eta.'*N2, xbsa,xbsa11,xbsa22),eta,-1,1));

Ms=1/Gs*double(int(int((N2.')*a11*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
Cs=1/Gs*double(int(int((N2.')*a22*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
Ks=1/Gs*double(int(int((N2.')*a33*N2, xbsa,xbsa11,xbsa22),eta,-1,1));
Plusor=diag(Mu);
if all(Plusor(:)<0)==1
    Mu=-Mu;
    Cu=-Cu;
    Ku=-Ku;
    CuS=-CuS;
    KuS=-KuS;
end
% sum()

% CSu=CuS.';
% % 
% % KSu2=double(   int(  int( N2.'*Fb*B_eta*Ea*B_xbsa_eta, xbsa,xbsa11,xbsa22),eta,-1,1));
% 
% KSu=KuS.';


% DeDe=[2*nodXY11(:,1)-1;nodXY22(:,1)];
de1=[nodXY11(:,1);nodXY22(:,1)];
dd=[[1:4]',de1];
de11=sortrows(dd,2);
de=de11(:,1);
De=zeros(8,1);
De(1:2:end,1)=2*de-1;
De(2:2:end,1)=2*de;

MU=Mu(de,de);
CU=Cu(de,de);
KU=Ku(de,de);
CUS=CuS(de,De);
KUS=KuS(de,De);

MS=Ms(De,De);
CS=Cs(De,De);
KS=Ks(De,De);
CSU=CSu(De,de);
KSU=KSu(De,de);

% CCC=Inv_Yakobi*B_xbsa_eta;
% aaa=Yakobi_det*ShapeFunM*rs.'*ShapeFunM.';
% MM=double(int(int(aaa, eta,-1,1),xbsa,xbsa11,xbsa22));
% MM=MM(de,de);

% 
% KK1=KK(de,:);
% 
% CSU=CUS.';
% KSU=KUS.';
% KS=Ks(De,De);


% Ku=ny*rs*double(int(int((Shape).'*c1*Shape, r,r3,r1),s,s3,s1));


% double(subs(sum(ShapeFunM(:,1),1),[xbsa,eta],[10,0.3]));


% esp=1.0e-5;
% for jj=1:4
%     for ii=1:4
%         Value=(eval(subs(ShapeFunM(ii,1),[r, s],[Coor_rs(1,jj) ,Coor_rs(2,jj) ])));
%         if  abs(Value -1)<esp
%             Shape(1,jj)=ShapeFunM(ii,1);
%         end
%     end
% end

% a_coefficient=(eval(subs(a1,[x y],[xx0,yy0])));
% 
% dy./dx
% 
% x0
% y0



% 
% 
% XX=(X1+X2+X3+X4)/4;
% YY=(Y1+Y2+Y3+Y4)/4;
% syms x y
% beta0=(m+1)*C_ref/(2*Lpml)*log(1/abs(R));
% alfa0=(m+1)*b_ele/(2*Lpml)*log(1/abs(R));
% 
% if abs(XX)<x0 && abs(YY)>y0
%     ns=1;
%     alfa_x=1;
%     beta_x=0;
%     beta_y=beta0*((y-y0)*ns/Lpml)^m;
%     alfa_y=1+alfa0*((y-y0)*ns/Lpml)^m;
% elseif abs(XX)>x0 && abs(YY)<y0
%     if XX>0;
%         ns=1;
%         alfa_x=1+alfa0*((x-x0)*ns/Lpml)^m;
%         beta_x=beta0*((x-x0)*ns/Lpml)^m;
%     else
%         ns=-1;
%         alfa_x=1+alfa0*((x+x0)*ns/Lpml)^m;
%         beta_x=beta0*((x+x0)*ns/Lpml)^m;
%     end
%     beta_y=0;
%     alfa_y=1;
% elseif abs(XX)>x0 && abs(YY)>y0
%     if XX>0;
%         ns=1;
%         alfa_x=1+alfa0*((x-x0)*ns/Lpml)^m;
%         beta_x=beta0*((x-x0)*ns/Lpml)^m;
%     else
%         ns=-1;
%         alfa_x=1+alfa0*((x+x0)*ns/Lpml)^m;
%         beta_x=beta0*((x+x0)*ns/Lpml)^m;
%     end
%     ns=1;
%     beta_y=beta0*((y-y0)*ns/Lpml)^m;
%     alfa_y=1+alfa0*((y-y0)*ns/Lpml)^m;
% end
% 
% % if abs(XX)<x0 && abs(YY)>y0
% %     ns=1;
% %     alfa_x=1;
% %     beta_x=0;
% %     beta_y=beta0*((y-y0)*ns/Lpml)^m;
% %     alfa_y=1+alfa0*((y-y0)*ns/Lpml)^m;
% % elseif abs(XX)>x0 && abs(YY)<y0
% %     if XX>0;
% %         ns=1;
% %         alfa_x=1+alfa0*((x-x0)*ns/Lpml)^m;
% %         beta_x=beta0*((x-x0)*ns/Lpml)^m;
% %     else
% %         ns=-1;
% %         alfa_x=1+alfa0*((x+x0)*ns/Lpml)^m;
% %         beta_x=beta0*((x+x0)*ns/Lpml)^m;
% %     end
% %     beta_y=0;
% %     alfa_y=1;
% % elseif abs(XX)>x0 && abs(YY)>y0
% % if XX>0;
% %     ns=1;
% %     alfa_x=1+alfa0*((x-x0)*ns/Lpml)^m;
% %     beta_x=beta0*((x-x0)*ns/Lpml)^m;
% % else
% %     ns=-1;
% %     alfa_x=1+alfa0*((x+x0)*ns/Lpml)^m;
% %     beta_x=beta0*((x+x0)*ns/Lpml)^m;
% % end
% % ns=1;
% % beta_y=beta0*((y-y0)*ns/Lpml)^m;
% % alfa_y=1+alfa0*((y-y0)*ns/Lpml)^m;
% % end
% 
% % PML 系数 alfa_y
% % alfa_x=1; alfa_y=1;
% a1=alfa_x*alfa_y;
% b1=alfa_x*beta_y+beta_x*alfa_y;
% C1=beta_x*beta_y;
% %  [SE_PML,ME_PML] = PML_M(n,L,Lp,Fre,cs,Gs,NOC_PML,rs)
% MA = zeros (4,4); MB = zeros (4,4); MC = zeros (4,4);
% A_alfay= zeros (4,4); A_alfax = zeros (4,4);
% A_betay = zeros (4,4);A_betax = zeros (4,4);
% 
% XNI(1, 1) = -0.57735026919;
% XNI(2, 1) = -0.57735026919;
% XNI(1, 2) =  0.57735026919;
% XNI(2, 2) = -0.57735026919;
% XNI(1, 3) =  0.57735026919;
% XNI(2, 3) =  0.57735026919;
% XNI(1, 4) = -0.57735026919;
% XNI(2, 4) =  0.57735026919;
% for p = 1 :4;
%     XI = XNI(1,p);
%     ETA = XNI(2,p);
%     TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4;
%     TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4;
%     TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4;
%     TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4;
%     DJ = TJ11 * TJ22 - TJ12 * TJ21;
%     %  --- Inversion of JACOBIAN
%     A(1, 1) = TJ22 / DJ;
%     A(1, 2) = -TJ12 / DJ;
%     A(2, 1) = -TJ21 / DJ;
%     A(2, 2) = TJ11 / DJ;
%     %  --- Local derivatives of shape function
%     g(1, 1) = -(1 - ETA) / 4;
%     g(2, 1) = -(1 - XI) / 4;
%     g(1, 2) = (1 - ETA) / 4;
%     g(2, 2) = -(1 + XI) / 4;
%     g(1, 3) = (1 + ETA) / 4;
%     g(2, 3) = (1 + XI) / 4;
%     g(1, 4) = -(1 + ETA) / 4;
%     g(2, 4) = (1 - XI) / 4;
%     %  --- B = A * g
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
%     Bx=B(1,:);     By=B(2,:);
%     %  --- BB = pose the B in the strainBB
%     % Element Stiffness Matrix  SE
%     ShapeFunM(1,1)=0.25*(1-XI)*(1-ETA);
%     ShapeFunM(1,2)=0.25*(1+XI)*(1-ETA);
%     ShapeFunM(1,3)=0.25*(1+XI)*(1+ETA);
%     ShapeFunM(1,4)=0.25*(1-XI)*(1+ETA);
%     xx0=ShapeFunM*[X1;X2;X3;X4];
%     yy0=ShapeFunM*[Y1;Y2;Y3;Y4];
%     
%     
%     a_coefficient=(eval(subs(a1,[x y],[xx0,yy0])));
%     b_coefficient=(eval(subs(b1,[x y],[xx0,yy0])));
%     c_coefficient=(eval(subs(C1,[x y],[xx0,yy0])));
%     
%     ay_coefficient=eval(subs(alfa_y,y,yy0));
%     ax_coefficient=eval(subs(alfa_x,x,xx0));
%     by_coefficient=eval(subs(beta_y,y,yy0));
%     bx_coefficient=eval(subs(beta_x,x,xx0));
%     
%     MA=MA+a_coefficient*ShapeFunM'*ShapeFunM*DJ;
%     MB=MB+b_coefficient*ShapeFunM'*ShapeFunM*DJ;
%     MC=MC+c_coefficient*ShapeFunM'*ShapeFunM*DJ;
%     
%     A_alfay=A_alfay+ay_coefficient*Bx.'*ShapeFunM*DJ;
%     A_alfax=A_alfax+ax_coefficient*By.'*ShapeFunM*DJ;
%     
%     A_betay=A_betay+by_coefficient*Bx.'*ShapeFunM*DJ;
%     A_betax=A_betax+bx_coefficient*By.'*ShapeFunM*DJ;
%     
% end
% A_alfay1=A_alfay.'*Gs;
% A_alfax1=A_alfax.'*Gs;
% A_betay1=A_betay.'*Gs;
% A_betax1=A_betax.'*Gs;
% 
% MA1=rs*MA;
% % MB1=rs*MB;
% MC1=rs*MC;
end

