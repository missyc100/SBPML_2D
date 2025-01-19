function [ UU,Dis] = newmark( MM,CC,KK,ntime,dt,NNs,Fix,Force_N,pressure)
N_sum=size(MM,1);
Feq=zeros(N_sum,ntime);
Feq(Force_N(2:end-1,1),:)=ones(size(Force_N(2:end-1,1),1),1)*pressure;
Feq(Force_N(1,1),:)=1/2*ones(size(Force_N(1,1),1),1)*pressure;
Feq(Force_N(end,1),:)=1/2*ones(size(Force_N(end,1),1),1)*pressure;

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
Kps_eq=KK+a0*MM+a1*CC;
Kps_eq(Fix,:)=0;
Kps_eq(:,Fix)=0;
% for jj=1:length(NB)
%     Kps_eq(NB(jj),NB(jj))=1;
% end
    Kps_eq(Fix,Fix)=eye(length(Fix));
UU=zeros(N_sum,ntime);
VV=zeros(N_sum,ntime);
AA=zeros(N_sum,ntime);
Feff = zeros(N_sum , ntime );
Keq=inv(Kps_eq);
% sparse(Kps_eq)
for j = 2 : ntime
    time=(j-1)*dt;
    fprintf('time %f **************\r\n',time);
    Feff(:,j) = Feq(:,j) + MM*(a0.*UU(:,j-1) + a2.*VV(:,j-1) + a3.*AA(:,j-1)) + CC*(a1.*UU(:,j-1) + a4.*VV(:,j-1) + a5.*AA(:,j-1));
    Feff(Fix,j)=0;
    UU(:,j) = Keq*Feff(:,j);
    AA(:,j) = a0.*(UU(:,j)-UU(:,j-1)) - a2.*VV(:,j-1) - a3.*AA(:,j-1);
    VV(:,j) = VV(:,j-1) + a6.*AA(:,j-1) + a7.*AA(:,j);
end
Dis=UU(1:2*NNs,:).';
% Stres_x=VV(NNs+1:2*NNs,:).';
% Stres_y=VV(2*NNs+1:3*NNs,:).';
end

