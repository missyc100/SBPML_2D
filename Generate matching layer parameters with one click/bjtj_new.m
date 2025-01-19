function [Uf]=bjtj_new(Fc,CoulpeK,plus_BL,Stress_free_Lxz,Stress_free_Lyz,Free_ux,...
    Absorb_matrices,IU_free_B,Fre,Stress_in_Byz,nxL,a_B1_inf,XofBottom,...
    Stress_In_Rxz,a_B2_inf,Sb_R,Plus_BRR,IU_free_R,mmd,nxR,Sb_L,cx,ARR,ARS)
nns=size(CoulpeK);
Feqr=zeros(nns(1),1); %Feqr=zeros(2*NB,1); 
% Feqr(Plus_BRR,1)=(-1i*Fre/cx*ARR+ARS)*IU_free_R+Sb_R*IU_free_R;

R_Num_node=size(a_B2_inf,1);
NOCR=zeros(R_Num_node-1,2);
NOCR(:,1)=1:R_Num_node-1; NOCR(:,2)=2:R_Num_node;
R_dy_mesh0=diff(a_B2_inf(:,2:3));
R_dy_mesh=sqrt(R_dy_mesh0(:,1).^2+R_dy_mesh0(:,2).^2);
Frz=zeros(R_Num_node,1); 
for ele=1:R_Num_node-1
    Fredom1=NOCR(ele, 1); Fredom2=NOCR(ele, 2);
    Shape_R=[1/3*R_dy_mesh(ele),1/6*R_dy_mesh(ele);1/6*R_dy_mesh(ele),1/3*R_dy_mesh(ele)];
    Frz(Fredom1:Fredom2,:)=Frz(Fredom1:Fredom2,:)+Shape_R*(nxR(ele,2)*Stress_In_Rxz(Fredom1:Fredom2,:));
end
Frz(end,1)=Frz(end,1)+Fc(1,1);
Frz=[Frz;Fc(2:end,1)];
Feqr(Plus_BRR,1)=Frz(:,1)+Sb_R*IU_free_R;

% % % FR=zeros(2*L_Num_node,1);
% % % FR(1:2:end)=Fx_R;FR(2:2:end)=Fy_R;
% % % 
Feql=zeros(nns(1),1);
L_dy_mesh0=diff(a_B1_inf(:,2:3));
L_dy_mesh=sqrt(L_dy_mesh0(:,1).^2+L_dy_mesh0(:,2).^2);
L_Num_node=size(a_B1_inf,1);
NOCL=zeros(L_Num_node-1,2);
NOCL(:,1)=1:L_Num_node-1; NOCL(:,2)=2:L_Num_node;
Fx_L=zeros(L_Num_node,1);
for ele=1:L_Num_node-1
    Fredom1=NOCL(ele, 1); Fredom2=NOCL(ele, 2);
    Shape_L=[1/3*L_dy_mesh(ele),1/6*L_dy_mesh(ele);1/6*L_dy_mesh(ele),1/3*L_dy_mesh(ele)];
    Fx_L(Fredom1:Fredom2,:)=Fx_L(Fredom1:Fredom2,:)+Shape_L*(-nxL(ele,2)*Stress_free_Lxz(Fredom1:Fredom2,:));
end

Feql(plus_BL,1)=Fx_L(:,1)+Sb_L*Free_ux;
% Feqr(plus_BR,1)=Sb_R*IU_free+ARS*IU_free+ARR*(1i*Fre*IU_free)./cx+BB*(Fre^2*IU_free)./((cx)^2);

% FF=FL;
% FS=-Sb_L*U_free_L(:,1);
% µ×²¿
Feqb=zeros(nns(1),1);
Num_node=size(XofBottom,1);
NOCbottom=zeros(Num_node-1,2);
NOCbottom(:,1)=1:Num_node-1; NOCbottom(:,2)=2:Num_node;
Fy_bottom=zeros(Num_node,1);
Bo_dx_mesh=diff(XofBottom(:,2));

for ele=1:Num_node-1
    Fredom1=NOCbottom(ele, 1); Fredom2=NOCbottom(ele, 2);
    Shape_bottom=[1/3*Bo_dx_mesh(ele),1/6*Bo_dx_mesh(ele);1/6*Bo_dx_mesh(ele),1/3*Bo_dx_mesh(ele)];
    Fy_bottom(Fredom1:Fredom2,:)=Fy_bottom(Fredom1:Fredom2,:)+Shape_bottom*(Stress_in_Byz(Fredom1:Fredom2,:));
end
Fby=Fy_bottom(1:end,1)+1i*Fre*Absorb_matrices*IU_free_B(1:end,1);
Feqb(XofBottom(1:end,1),1)=+Fby(1:end,1);
Feq=Feql+Feqb+Feqr;
Uf=sparse(CoulpeK)\Feq(:,1);
end

