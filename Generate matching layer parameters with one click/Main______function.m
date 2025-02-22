clc
clear all
%%  全局变量
global  NNs NEs  NOC_Solid  XofNs  ...
    Gs rs nxL nxR  a_B1_inf a_B2_inf
%%  内域几何物理信息
NNs=561;
NEs=512;
NOC_Solid=dlmread('element.txt',',',[1,1,NEs,4]);
XofNs=dlmread('node.txt',',',[1,1,NNs,2]);
esp=1.0e-8;
% 求动力学矩阵
cs=2;rs=20;Gs=cs^2*rs; 
Lpml=input('下部匹配层厚度===  ');
R=1.0e-9;
b_ele=0.02; 

dt=0.002;ntime=5001; t0=0.8;fp=1.5;A_max=1;
t=(0:dt:(ntime-1)*dt);Tp=(ntime-1)*dt;
F0=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
F0(abs(F0)<esp)=0;
omiga=0:2*pi/(Tp):2*pi/(dt);
F0_fre =fft(F0);
plot(omiga./2./pi,abs(F0_fre))
Nt=2001;

m=2;
x0=0.4;y0=0.4;
[MM1,CC1,KK1] = Assembly(NOC_Solid,XofNs,Gs,rs,Lpml,R,b_ele,m,x0,y0);%% 边界处理
[ UU,Dis,Stres_x,Stres_y ] = newmark( MM1,CC1,KK1,F0,Nt,dt,NNs);
%% 结果
CoorxL=-0.12; CoorxR=0.12;Coory=0;
ro_imposeXL=find(abs(XofNs(:,1)-CoorxL)<1.0e3*esp);  
ro_imposeXR=find(abs(XofNs(:,1)-CoorxR)<1.0e3*esp);  
ro_imposeX=[ro_imposeXL;ro_imposeXR];
ro_imposeY=find(abs(XofNs(ro_imposeX,2)-Coory)<1.0e3*esp);  
xxx=ro_imposeX(ro_imposeY);
Uy_f=Dis(:,xxx);
result=[t(1,1:Nt).',Uy_f(:,1),t(1,1:Nt).',Uy_f(:,2)];


% % [N_B1,b]=find(abs(XofNs(:,1)-0)<esp);
% N_B1= importdata('BLL_node.txt');
% N_B1= N_B1(find( ~ isnan(N_B1)));
% Xof_B1=XofNs(N_B1,:);
% % % N_B1=[21:21:861]';
% % % Xof_B1=XofNs(N_B1,:);
% B1_inf=[N_B1,Xof_B1];
% a_B1_inf=sortrows(B1_inf,3);
% ele_vec=diff(a_B1_inf(:,2:3));
% unit_vec=zeros(length(ele_vec),2);
% for i=1:length(ele_vec)
%     unit_vec(i,:)=ele_vec(i,:)./ norm(ele_vec(i,:));
% end
% nx1=unit_vec*[1,0]'; ny1=unit_vec*[0,1]';
% nxL=[nx1,ny1];
% plus_BL=a_B1_inf(:,1);
% BNL=size(a_B1_inf,1);
% BEL=size(nxL,1);
% NOC_BL=[1:BNL-1;2:BNL]';
% % 单元离散矩阵
% [LArr0,LArs0,LAss0,LMM0] = Integrate_Solid(BNL,NOC_BL,BEL,a_B1_inf,nxL);
% Zl=zeros(BNL,BNL);
% %%  Right 边界部分
% N_B2= importdata('BRR_node.txt');
% N_B2= N_B2(find( ~ isnan(N_B2)));
% Xof_B2=XofNs(N_B2,:);
% B2_inf=[N_B2,Xof_B2];
% a_B2_inf=sortrows(B2_inf,3);
% depth=a_B2_inf(end,3);
% ele_vec2=diff(a_B2_inf(:,2:3));
% unit_vec2=zeros(length(ele_vec2),2);
% for i=1:length(ele_vec2)
%     unit_vec2(i,:)=ele_vec2(i,:)./ norm(ele_vec2(i,:));
% end
% nx2=unit_vec2*[1,0]'; ny2=unit_vec2*[0,1]';
% nxR=[nx2,ny2];
% BNR=size(a_B2_inf,1);
% plus_BR=a_B2_inf(:,1);
% BER=size(nxR,1);
% NOC_BR=[1:BNR-1;2:BNR]';
% [RArr0,RArs0,RAss0,RMM0] = Integrate_Solid(BNR,NOC_BR,BER,a_B2_inf,nxR); 
% Zr=zeros(BNR,BNR);
% %% 底部半空间物理信息
% GH =dlmread('Mate_half space.txt',',',[1,0,1,0]);
% denH =dlmread('Mate_half space.txt',',',[1,1,1,1]);
% CsH= sqrt(real(GH)/denH);
% [Ah0,Bh0,Ch0] = Augmented_matrix(CsH,GH);
%% 荷载
% dt=0.01;ntime=2001; t0=0.6;fp=2;A_max=1;
% t=(0:dt:(ntime-1)*dt);Tp=(ntime-1)*dt;
% loading=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
% loading(abs(loading)<esp)=0;
% omiga=0:2*pi/(Tp):2*pi/(dt);
% F0 =fft(loading);
% plot(omiga./2./pi,abs(F0))
% xxx=[t',loading'];
% Fre_max=fp*5*2*pi;
% [a,b]=find(omiga<Fre_max);
% max_nfre=max(b); %需要计算
% [CoorXY,pressure] = enforce_lump(F0);
% %% 存贮向量
% UU=zeros(NNs,ntime);
% Loading=zeros(NNs,ntime);
% Loading(CoorXY,:)=2*pressure;
% 
% % 注意每次循环需要注意， 即输出与输入都是什么值，还有下次循环还是不是初值
% for jj=2:max_nfre
%     Fre=omiga(jj);
%     % 计算边界
%     [LARR,LARS_RST,LASS_M,LARS] = Dry_Co(LArr0,LArs0,LAss0,LMM0,Fre);
%     % 基岩
%     if Fre~=0
%         Ah=Ah0/(1i*Fre);
%         Ch=(1i*Fre)*Ch0;
%         [AAL,ACL]=plus_Reig(LARR,LASS_M,Ah,Ch,BNL);
%     else
%         AAL=LARR;
%         ACL=LASS_M;
%         [LARS_RST,LARS]=plus_Reig0(LARS_RST,Bh0,LArs0,BNL);
%     end
%     ABL=1i*LARS_RST;
%     MBL=[-ACL,Zl;Zl,AAL];
%     KBL=[Zl,-ACL;-ACL,-ABL];
%     AeL=MBL\KBL;
%     [FeiL,ladL]=eig(AeL);
%     [LaL,FbL,jj4L]=solve_eig_L(FeiL,ladL,BNL,esp);
%     Sb_L=-(AAL*FbL*(1i)*LaL/(FbL)-LARS);
%     
%     % Right
%     [RARR,RARS_RST,RASS_M,RARS] = Dry_Co( RArr0,RArs0,RAss0,RMM0,Fre);
%     if Fre~=0
%         [AAR,ACR]=plus_Reig(RARR,RASS_M,Ah,Ch,BNR);
%     else
%         AAR=RARR;
%         ACR=RASS_M;
%         [RARS_RST,RARS]=plus_Reig0(RARS_RST,Bh0,RArs0,BNR);
%     end
%     ABR=1i*RARS_RST;
%     MBR=[-ACR,Zr;Zr,AAR];
%     KBR=[Zr,-ACR;-ACR,-ABR];
%     AeR=KBR\MBR;
%     [FeiR,ladR]=eig(AeR);
%     [LaR,FbR,jj4R]=solve_eig_R(FeiR,1./ladR,BNR,esp);
%     Sb_R=(AAR*FbR*(1i)*LaR/(FbR)-RARS);
%     
%     % 耦合
%     [CoulpeK]=K_boundary1(MM,KK,Sb_R,plus_BR,Sb_L,plus_BL,Fre);
%     UU(:,jj)=sparse(CoulpeK)\Loading(:,jj);
%     %       Ub(:,jj)=Sb_R\Loading(:,jj);
% end

% ro_imposeY=find(abs(XofNs(:,2)-1)<1.0e2*esp);
% ro_imposeX=find(abs(XofNs(ro_imposeY,1)-0)<1.0e2*esp);
% CoorXY=ro_imposeY(ro_imposeX);
% Uy_f=UU(CoorXY,:).';
% Uy_t2=(ifft(Uy_f,'symmetric'));
% plot(t',(Uy_t2),'--');
% hold on
% Reference solution
% co=xx_coor/ele_Y+1;


Ut2=UU(1:NNs,:).';
fid=fopen('U3y.postr.res','wt');
achar='GiD Post Results File 1.0';
bchar='result "Pressure" "2D Scalar Waveguide"   ';
Time=[0.2:0.01:2.4];%:0.02:10;
Co=round(Time/dt)+1;
Node=[(1:NNs)',XofNs];
ele=[(1:NEs)',NOC_Solid];
fprintf(fid,'%s\n',achar);


for jj=1:size(Time,2)
    ux_everynode=Ut2(Co(1,jj),1:end)';
    No_value=[(1:NNs)',ux_everynode];
    Tim=num2str(Time(1,jj));
    fprintf(fid,'%s',bchar);
    fprintf(fid,'%s',Tim);
    fprintf(fid,'%s','s');
    fprintf(fid,'%s\n','     Scalar OnNodes');
    fprintf(fid,'%s\n','values');
    num=size(No_value,1);
    for ii=1:num
        fprintf(fid,'%6.0d',No_value(ii,1));
        fprintf(fid,'%s','     ');
        fprintf(fid,'%d\n',No_value(ii,2));
    end
    fprintf(fid,'%s\n','end values');
end


fclose(fid);   



