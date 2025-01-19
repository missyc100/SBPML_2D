function [U1x,U1y,cx,A1x,A1y]...
    =OneDimensioF0n_FE_AB(Y_s,angle,PorSV,dt,ntime,F0)
% Coupling condition
% Bedrock half-space problem
% lamdH(1,1)= dlmread('material_constant_Half.txt',',',[1,0,1,0]);
% GH(1,1)= dlmread('material_constant_Half.txt',',',[1,1,1,1]);
% rsH(1,1)= dlmread('material_constant_Half.txt',',',[1,2,1,2]);
lamdH(1,1)=270000000;% Lamd constant
GH(1,1) = 180000000;% Shear modulus
rsH(1,1) = 2000;% density

cpH=sqrt((lamdH+2*GH)/rsH);% P-wave velocity
csH=sqrt((GH)/rsH);% SV-wave velocity
if PorSV==1
    cx=cpH(end,1)/sind(angle);% The incident wave is the apparent wave velocity of P wave
else PorSV==2
    cx=csH(end,1)/sind(angle);% The incident wave is the apparent wave velocity of SV wave
end
npx=cpH(end,1)/cx;
npy=sqrt(1-npx^2);
nsx=csH(end,1)/cx;
nsy=sqrt(1-nsx^2);
[S] = ABC(lamdH,GH,rsH(end,1),npx,npy,nsx,nsy,cx,PorSV);
% % Dry Soil
depth1=max(Y_s);
NN_S=length(Y_s);
NE_S=NN_S-1;
XofNS=zeros(NN_S,2);XofNS(:,1)=1:NN_S;XofNS(:,2)=Y_s;
NOCS=zeros(NE_S,2);NOCS(:,1)=1:NE_S;NOCS(:,2)=1:NN_S-1;NOCS(:,3)=2:NN_S;
lamds=ones(NE_S,1);Gs=ones(NE_S,1);rs=ones(NE_S,1);
Mid_ele=Y_s(1:end-1,:)+diff(Y_s)./2;
Drysoil_layer=1200;
for ji=1:NE_S
        lamds(ji,1)= dlmread('material_constant_Solid.txt',',',[ji,0,ji,0]);
        Gs(ji,1)= dlmread('material_constant_Solid.txt',',',[ji,1,ji,1]);
        rs(ji,1)= dlmread('material_constant_Solid.txt',',',[ji,2,ji,2]);
%     lamds(ji,1)=4.8*10^8;
%     Gs(ji,1)=3.2*10^8;
%     rs(ji,1)=2000;
%     lamds(ji,1)=9.991321200000000*10^9;
%     Gs(ji,1)=1.000188000000000*10^10;
%     rs(ji,1)=2800;

%     z = max(Y_s)-(Y_s(ji,1)+Y_s(ji,1)+1)/2;
%     z = max(Y_s)-Y_s(ji+1,1);
%     Vh = 400;%最大剪切波速
%     V0 = 200;%最小剪切波速
%     n = 0.5;
%     b = (V0/Vh)^n;
%     H = max(Y_s)-min(Y_s);
%     Vs = Vh*(b+(1-b)*z/H)^n;
%     rs(ji,1)=2000;
%     v = 0.3;
%     Gs(ji,1)=Vs^2*rs(ji,1);
%     lamds(ji,1)=2*Gs(ji,1)*v/(1-2*v);
end
 

[AMs,ACs,AKs]=integrate_S(XofNS,NOCS,cx,lamds,Gs,rs);

N_sum = 2*NN_S;
Mps=zeros(N_sum,N_sum);
Cps=zeros(N_sum,N_sum);
Kps=zeros(N_sum,N_sum);
%%%%形成内域总刚阵和总质量阵
Mps(1:2*NN_S,1:2*NN_S)=AMs;

Cps(1:2*NN_S,1:2*NN_S)=ACs;

Kps(1:2*NN_S,1:2*NN_S)=AKs;
% 三者耦合


N_sum=length(Cps);
% %%%%%%%荷载%%%%%%%%%%%%%%%(换荷载时要改）
% dt = 0.001;
% ntime =4001;
% t = dt.*(0:1:ntime-1)';
% A =0.01;
% T =0.5;
% x =(t)./T;
% u = A.*16.*(x.^3.*heaviside(x)-4.*(x-0.25).^3.*heaviside(x-0.25)+6.*(x-0.5).^3.*heaviside(x-0.5)-4.*(x-0.75).^3.*heaviside(x-0.75)+(x-1).^3.*heaviside(x-1));
% F0 = A.*16.*3./T.*(x.^2.*heaviside(x)-4.*(x-0.25).^2.*heaviside(x-0.25)+6.*(x-0.5).^2.*heaviside(x-0.5)-4.*(x-0.75).^2.*heaviside(x-0.75)+(x-1).^2.*heaviside(x-1));
% plot(t,u)
% Ear=xlsread('History_Earthquake.xlsx');
%% 入射波时程——————需要换
% esp=1.0e-8;
% t0=1;fp=2;A_max=0.001;
% t=(0:dt:(ntime-1)*dt);Tp=(ntime-1)*dt;
% loading=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
% loading(abs(loading)<esp)=0;  % 位移
% % plot(t,loading)
% % incidence=fft(loading);
% loading1=(-4*A_max*pi^2*fp^2*(t-t0)-2*A_max*pi^2*fp^2*(t-t0).*(1-2*(pi*fp*(t-t0)).^2)).*exp(-(pi*fp*(t-t0)).^2);
% loading1(abs(loading1)<esp)=0; % 速度
% % omiga=(0:2*pi/(Tp):2*pi/(dt));
% % plot(omiga./(2*pi),abs(incidence))
% F0 =loading1';%fft(loading);  
% Fre_max=fp*10*2*pi;
%% %%%%%%%%%%%%%%%%%%%%%%%%%
% [a,b]=find(omiga<Fre_max);
% max_nfre=max(b); %需要计算

% plot(omiga./(2*pi), abs(F0))
% ntime=ntime/space_n;
%%%%%%%%%%%定义U-u 方程中的两个耦合方程%%
% Newmark coeffient开始计算
bata=1/4;
gama=1/2;
a0=1/(bata*dt^2);
a1 = gama/(bata*dt);
a2 = 1/(bata*dt);
a3 = 1/(2*bata) - 1;
a4 = gama/bata - 1;
a5 = dt/2*(gama/bata - 2);
a6 = dt * (1 - gama);
a7 = gama*dt;
% Boundary condition
if PorSV == 1
    F2 = ( S * [npx(1);-npy(1)] + [2*GH(1)*npx(1)*npy(1);-lamdH(1)-2*GH(1)*npy(1)^2] ./ cpH(1) ) * F0'; % the load of artificial boundary
elseif PorSV == 2
    F2 = ( S * [nsy(1);nsx(1)] +[-(nsx(1)^2-nsy(1)^2)*GH(1); 2*GH(1)*nsx(1)*nsy(1)]./csH(1).*(1) ) * F0'; % the load of artificial boundary
end
Cps(N_sum-1:N_sum,N_sum-1:N_sum)=Cps(N_sum-1:N_sum,N_sum-1:N_sum)+S;% the equivalent stiffness matrice
Kps_eq=Kps+a0*Mps+a1*Cps;
% 干土与饱和间透水
%%%%%%%%%%%%%Newmark_bata法求解计算%%%%%
UU=zeros(N_sum,ntime);
VV=zeros(N_sum,ntime);
AA=zeros(N_sum,ntime);
F = zeros(N_sum , ntime );
Feff = zeros(N_sum , ntime );
F(N_sum-1:N_sum,:)=  F(N_sum-1:N_sum,:)+F2(1:2,:);
for j = 2 : ntime
    %        time=(j-1)*dt;
    %     fprintf('time %f **************\r\n',time);
    Feff(:,j) = F(:,j) + Mps*(a0.*UU(:,j-1) + a2.*VV(:,j-1) + a3.*AA(:,j-1)) + Cps*(a1.*UU(:,j-1) + a4.*VV(:,j-1) + a5.*AA(:,j-1));
    UU(:,j) = sparse(Kps_eq)\Feff(:,j);
    AA(:,j) = a0.*(UU(:,j)-UU(:,j-1)) - a2.*VV(:,j-1) - a3.*AA(:,j-1);
    VV(:,j) = VV(:,j-1) + a6.*AA(:,j-1) + a7.*AA(:,j);
end
%% 后处理 output displacement

% 将公用节点自由度位置补上
Co_UU =UU;
Co_VV =VV;
Co_AA = AA;
% 将不同介质的动力响应分离：
% 1—dry soil;
u1x = Co_UU(1:2:2*NN_S,:);    u1y = Co_UU(2:2:2*NN_S,:);
v1x = Co_VV(1:2:2*NN_S,:);    v1y = Co_VV(2:2:2*NN_S,:);
a1x = Co_AA(1:2:2*NN_S,:);    a1y = Co_AA(2:2:2*NN_S,:);
U1x=u1x'; U1y=u1y';
A1x=a1x'; A1y=a1y';
% % Dry soil应力计算
% U_solid= Co_UU(1:2*NN_S,:);V_solid= Co_VV(1:2*NN_S,:);A_solid= Co_AA(1:2*NN_S,:);
% Sigma_S=zeros(2*NN_S,ntime);
% for ele=1:NE_S
%     Fredom1=NOCS(ele, 2);
%     Fredom2=NOCS(ele, 3);
%     xy1=2*Fredom1-1:2*Fredom1;
%     xy2=2*Fredom2-1:2*Fredom2;
%     xy=[xy1,xy2];
%     Sigma_S(xy,:)=MMs(ele,XofNS,NOCS,cx,lamds,Gs,rs)*A_solid(xy,:)+CCs(ele,XofNS,NOCS,cx,lamds,Gs,rs)*V_solid(xy,:)+KKs(ele,XofNS,NOCS,cx,lamds,Gs,rs)*U_solid(xy,:);
% end
% Sigma_S(1:end-2,:)=-Sigma_S(1:end-2,:);
% SigmaS_xy=Sigma_S(1:2:end,:);SigmaS_y=Sigma_S(2:2:end,:);
% lamds=([lamds(1);lamds]+[lamds;lamds(end)])./2;
% Gs=([Gs(1);Gs]+[Gs;Gs(end)])./2;
% SigmaS_x=(-diag(lamds+2*Gs)*V_solid(1:2:end,:)/cx+diag(lamds./(lamds+2*Gs)) *(SigmaS_y+diag(lamds)*V_solid(1:2:end,:)/cx));
% 
% 
% %% 快速Fourier Transformation
% Free_u1x=(fft(u1x')).';   Free_u1y=(fft(u1y')).';
% Fre_uxy=zeros(2*NN_S,ntime);
% Fre_uxy(1:2:end,:)=Free_u1x;  Fre_uxy(2:2:end,:)=Free_u1y;
% Free_SigmaS_x=(fft(SigmaS_x')).';   Free_SigmaS_y=(fft(SigmaS_y')).';   Free_SigmaS_xy=(fft(SigmaS_xy')).';
end

