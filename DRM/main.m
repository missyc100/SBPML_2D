close all;clear all;
clc
format long e
if exist('output1','dir')~=0
    rmdir('output1', 's')
end
if exist('output2','dir')~=0
    rmdir('output2', 's')
end
if exist('output4','dir')~=0
    rmdir('output4', 's')
end
pause(5)
node_all = readmatrix('input/node_all.txt');% abaqus model node
element_all = readmatrix('input/element_all.txt');% abaqus model elements
cross= importdata('input/cross.txt'); % Load the node number that the SBPML shares with the CPE4 elements
cross= cross(find( ~ isnan(cross))); 
cross = sort(cross);
DRM_ele=importdata('input/DRM_ele.txt');
DRM_ele= DRM_ele(find( ~ isnan(DRM_ele)));
DRM_ele = sort(DRM_ele);
DRM_ele = element_ex(element_all,DRM_ele);
node_e = cross;
node_b = setdiff(unique(DRM_ele(:,2:5)),node_e);
node_e = node_ex(node_all,node_e);
node_b = node_ex(node_all,node_b);

% % ricker
% dt=0.0002;
% T = 4;
% ntime=T/dt+1; A_max=1e-3;fp=3;t0=0.5;
% esp=1.0e-36;
% t=(0:dt:(ntime-1)*dt);Tp=(ntime-1)*dt;
% loading=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
% loading(abs(loading)<esp)=0;
% pp = loading';
% incidence=fft(loading);
% loading1=(-4*A_max*pi^2*fp^2*(t-t0)-2*A_max*pi^2*fp^2*(t-t0).*(1-2*(pi*fp*(t-t0)).^2)).*exp(-(pi*fp*(t-t0)).^2);
% loading1(abs(loading1)<esp)=0;  %% 

%地震动
% AA = readmatrix('G:\abaqus\paper1\DRM\app2\输入\northbridge_v.xlsx');
% AA(:,2) = AA(:,2)/200;%
% dt = AA(2,1)-AA(1,1);
% ntime = length(AA);
% A_max = max(AA(:,2));
% esp=1.0e-30;
% t = AA(:,1)';
% loading1 = (AA(:,2))';

% %dirac
dt=0.0002;
T = 3.0;
ntime=T/dt+1; 
ntime = fix(ntime);A=15;sigma=0.4/A;t0=0.3;A_max = 0.001;
esp=1.0e-10;
t=(0:dt:(ntime-1)*dt);Tp=(ntime-1)*dt;
loading=((1/sqrt(2*pi)./sigma).*exp(-(t-t0).^2/2/sigma^2))*A_max/A;
loading(abs(loading)<esp)=0;
pp = loading';
incidence=fft(loading);
loading1=(-sqrt(2).*(t-t0).*exp(-(t-t0).^2/(2*sigma^2))/(2*sigma^3*sqrt(pi)))*A_max/A;
loading1(abs(loading1)<esp)=0;  %% Incident wave velocity time history

Fei1=0;% Angle of incidence
Fei2=0;% The conversion Angle is 0 in the 2D case
YY=sortrows(unique([node_e(:,3);node_b(:,3)]),-1);% The height of nodes in DRM elements
% YY_add = min(YY):0.2:max(YY);
% YY_add = YY_add';
% YY = sortrows(unique([YY;YY_add]),-1);
%% Site reflection calculation
PorSV = 1;% 1 is a P wave and 2 is a SV wave
[U1x,U1y,cx2,A1x,A1y]=OneDimension_FE_AB(YY,Fei1,PorSV,dt,ntime,loading1.');
U1y(abs(U1y)<esp)=0;
U1x(abs(U1x)<esp)=0;
A1x(abs(A1x)<esp)=0;
A1y(abs(A1y)<esp)=0;
U1y = -U1y;
A1y = -A1y;

x_start = min(unique([node_b(:,2);node_e(:,2)]));% Displacement reference column

%% Displacement calculation of node_b nodes
for i1 = 1:length(node_b(:,1))
    ndt_b(i1) = (node_b(i1,2)-x_start)./cx2;
    [ro,co]=find(abs(YY(:,1)-node_b(i1,3))<esp);
    if ndt_b(i1) == 0
        Free_bx(i1,:) = (U1x(:,ro))';
        Free_by(i1,:) = (U1y(:,ro))';
        Free_bax(i1,:) = (A1x(:,ro))';
        Free_bay(i1,:) = (A1y(:,ro))';
    else
        tt = [];
        tt = t+ndt_b(i1);
        tt = [0,tt];
        dis_bx = [];
        dis_bx = (U1x(:,ro))';
        dis_bx = [0,dis_bx];
        Free_bx(i1,:) = interp1(tt,dis_bx,t);
        dis_by = [];
        dis_by = (U1y(:,ro))';
        dis_by = [0,dis_by];
        Free_by(i1,:) = interp1(tt,dis_by,t);
        dis_bax = [];
        dis_bax = (A1x(:,ro))';
        dis_bax = [0,dis_bax];
        Free_bax(i1,:) = interp1(tt,dis_bax,t);
        dis_bay = [];
        dis_bay = (A1y(:,ro))';
        dis_bay = [0,dis_bay];
        Free_bay(i1,:) = interp1(tt,dis_bay,t);
    end
end

%% node_e Node displacement calculation
for i1 = 1:length(node_e(:,1))
    ndt_e(i1) = (node_e(i1,2)-x_start)./cx2;
    [ro,co]=find(abs(YY(:,1)-node_e(i1,3))<esp);
    if ndt_e(i1) == 0
        Free_ex(i1,:) = (U1x(:,ro))';
        Free_ey(i1,:) = (U1y(:,ro))';
        Free_eax(i1,:) = (A1x(:,ro))';
        Free_eay(i1,:) = (A1y(:,ro))';
    else
        tt = [];
        tt = t+ndt_e(i1);
        tt = [0,tt];
        dis_ex = [];
        dis_ex = (U1x(:,ro))';
        dis_ex = [0,dis_ex];
        Free_ex(i1,:) = interp1(tt,dis_ex,t);
        dis_ey = [];
        dis_ey = (U1y(:,ro))';
        dis_ey = [0,dis_ey];
        Free_ey(i1,:) = interp1(tt,dis_ey,t);
        dis_eax = [];
        dis_eax = (A1x(:,ro))';
        dis_eax = [0,dis_eax];
        Free_eax(i1,:) = interp1(tt,dis_eax,t);
        dis_eay = [];
        dis_eay = (A1y(:,ro))';
        dis_eay = [0,dis_eay];
        Free_eay(i1,:) = interp1(tt,dis_eay,t);
    end
end

%% Form a node collection inp file
ampx = 1;ampy = 1;
ret = nodeset_inp(node_b(:,1),node_e(:,1),ampx,ampy);
%% The be stiffness matrix is formed
K_FileName = 'MTX\Job-KKMM_STIF2.mtx';
M_FileName = 'MTX\Job-KKMM_MASS3.mtx';
[KKbe,KKeb,MMbe,MMeb]=GENERATE_KK(node_b(:,1),node_e(:,1),K_FileName,M_FileName);
%% Time history inp -- include
dir = 'output2';
mkdir(dir)
node_index = [node_b(:,1);node_e(:,1)];
for i = 1:size(node_index,1)
    message1{1,1} = '*include';
    message1{1,2} = ['input=case1\output4\outputP',char(string(node_index(i,1))),'x.txt'];
    message2{1,1} = '*include';
    message2{1,2} = ['input=case1\output4\outputP',char(string(node_index(i,1))),'y.txt'];

    writecell(message1,[dir,'\include1.txt'],'WriteMode','append')
    writecell(message2,[dir,'\include1.txt'],'WriteMode','append')

end
%% Form a node load time history inp file
ret = node_timehistory_inp(Free_bx,Free_by,Free_ex,...
    Free_ey,node_b(:,1),node_e(:,1),KKbe,KKeb,dt,MMbe,MMeb,Free_bax,...
    Free_bay,Free_eax,Free_eay);







