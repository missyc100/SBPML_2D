function [uz,Free_ux,Free_Tx,Free_Ty,Ty,Ah,incidence]=OneDimension_FE_AB(Y_s,angle,dt,ntime,t0,fp,A_max)
% 耦合情况
% 基岩半空间问题
GH(1,1)= dlmread('Mate_half space.txt',',',[1,0,1,0]);
rsH(1,1)= dlmread('Mate_half space.txt',',',[1,1,1,1]);
cs=sqrt(GH./rsH);
cx=cs(end,1)/sind(angle);
Ah = rsH*cs*cosd(angle);
% %  Dry Soil
NN_S=length(Y_s);
NE_S=NN_S-1;
NOCS=zeros(NE_S,3);NOCS(:,1)=1:NE_S;NOCS(:,2)=1:NN_S-1;NOCS(:,3)=2:NN_S;

[MM,CC,KK]=integrate_S_1D(NN_S,NOCS(:,2:3),NE_S,Y_s,Ah,cx);


N_sum=length(CC);
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

esp=1.0e-9;
t=(0:dt:(ntime-1)*dt);Tp=(ntime-1)*dt;
loading=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
loading(abs(loading)<esp)=0;
% plot(t,loading)
incidence=fft(loading);
loading1=(-4*A_max*pi^2*fp^2*(t-t0)-2*A_max*pi^2*fp^2*(t-t0).*(1-2*(pi*fp*(t-t0)).^2)).*exp(-(pi*fp*(t-t0)).^2);
loading1(abs(loading1)<esp)=0;
omiga=(0:2*pi/(Tp):2*pi/(dt));
% plot(omiga./(2*pi),abs(incidence))
F0 =loading1';%fft(loading);  

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
F2 = 2*Ah* F0'; % the load of artificial boundary
Kps_eq=KK+a0*MM+a1*CC;
% 干土与饱和间透水
%%%%%%%%%%%%%Newmark_bata法求解计算%%%%%
UU=zeros(N_sum,ntime);
VV=zeros(N_sum,ntime);
AA=zeros(N_sum,ntime);
F = zeros(N_sum , ntime );
Feff = zeros(N_sum , ntime );
F(N_sum,:)=  F(N_sum,:)+F2;
for j = 2 : ntime
    %        time=(j-1)*dt;
    %     fprintf('time %f **************\r\n',time);
    Feff(:,j) = F(:,j) + MM*(a0.*UU(:,j-1) + a2.*VV(:,j-1) + a3.*AA(:,j-1)) + CC*(a1.*UU(:,j-1) + a4.*VV(:,j-1) + a5.*AA(:,j-1));
    UU(:,j) = sparse(Kps_eq)\Feff(:,j);
    AA(:,j) = a0.*(UU(:,j)-UU(:,j-1)) - a2.*VV(:,j-1) - a3.*AA(:,j-1);
    VV(:,j) = VV(:,j-1) + a6.*AA(:,j-1) + a7.*AA(:,j);
end
%% 后处理 output displacement
% 将公用节点自由度位置补上
% 将不同介质的动力响应分离：
% 1—dry soil;
Ty=zeros(NN_S,ntime);
for ele=1:NE_S
    Fredom1=NOCS(ele, 2);
    Fredom2=NOCS(ele, 3);
    xy1=Fredom1;
    xy2=Fredom2;
    xy=[xy1,xy2];
    [EA11,EA22,EMM] = Ele_matrix_1D(ele,Y_s);
    Ty(xy,:)=(EMM-EA11./cx^2)*AA(xy,:)+EA22*UU(xy,:);
end
Ty(1:end-1,:)=-Ty(1:end-1,:);
Tx=-GH/cx*VV;
%% 快速Fourier Transformation
Free_ux=(fft(UU')).'; 
Free_Tx=(fft(Tx')).';   Free_Ty=(fft(Ty')).';  
uz=UU';
end

