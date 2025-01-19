function [Uf]=bjtj_old(CoulpeK,Sb_R,BB,ARR,ARS,plus_BR,cx,Fre,IU_free,Stress_free_Lxz,Stress_free_Lyz,a_B2_inf,nxR)
nns=size(CoulpeK);
Feqr=zeros(nns(1),1); %Feqr=zeros(2*NB,1); %Feqb=zeros(2*NNs,1);
% L_Num_node=size(a_B1_inf,1);
R_Num_node=size(a_B2_inf,1);
% NOCL=zeros(L_Num_node-1,2);
% NOCL(:,1)=1:L_Num_node-1; NOCL(:,2)=2:L_Num_node;
NOCR=zeros(R_Num_node-1,2);
NOCR(:,1)=1:R_Num_node-1; NOCR(:,2)=2:R_Num_node;
% L_dy_mesh0=diff(a_B1_inf(:,2:3));
R_dy_mesh0=diff(a_B2_inf(:,2:3));
R_dy_mesh=sqrt(R_dy_mesh0(:,1).^2+R_dy_mesh0(:,2).^2);
% L_dy_mesh=sqrt(L_dy_mesh0(:,1).^2+L_dy_mesh0(:,2).^2);
% %
Fz=zeros(R_Num_node,1); 
for ele=1:R_Num_node-1
    Fredom1=NOCR(ele, 1); Fredom2=NOCR(ele, 2);
    Shape_R=[1/3*R_dy_mesh(ele),1/6*R_dy_mesh(ele);1/6*R_dy_mesh(ele),1/3*R_dy_mesh(ele)];
    Fz(Fredom1:Fredom2,:)=Fz(Fredom1:Fredom2,:)+Shape_R*(nxR(ele,2)*Stress_free_Lxz(Fredom1:Fredom2,:)-nxR(ele,1)*Stress_free_Lyz(Fredom1:Fredom2,:));
end
% % % FR=zeros(2*L_Num_node,1);
% % % FR(1:2:end)=Fx_R;FR(2:2:end)=Fy_R;
% % 
% Fx_L=zeros(L_Num_node,1); Fy_L=zeros(L_Num_node,1);
% for ele=1:L_Num_node-1
%     Fredom1=NOCL(ele, 1); Fredom2=NOCL(ele, 2);
%     Shape_L=[1/3*L_dy_mesh(ele),1/6*L_dy_mesh(ele);1/6*L_dy_mesh(ele),1/3*L_dy_mesh(ele)];
%     Fx_L(Fredom1:Fredom2,:)=Fx_L(Fredom1:Fredom2,:)+Shape_L*(-nxL(ele,2)*Stress_free_Lx(Fredom1:Fredom2,:)+nxL(ele,1)*Stress_free_Lxy(Fredom1:Fredom2,:));
%     Fy_L(Fredom1:Fredom2,:)=Fy_L(Fredom1:Fredom2,:)+Shape_L*(nxL(ele,1)*Stress_free_Ly(Fredom1:Fredom2,:)-nxL(ele,2)*Stress_free_Lxy(Fredom1:Fredom2,:));
% end
% FL=zeros(2*L_Num_node,1);
% FL(1:2:end)=Fx_L;FL(2:2:end)=Fy_L;
% Feql(plus_BLL,1)=FL(:,1)+Sb_L*U_free_L(:,1);

Feqr(plus_BR,1)=Fz(:,1)+Sb_R*IU_free(:,1);
% Feqr(plus_BR,1)=Sb_R*IU_free+ARS*IU_free+ARR*(1i*Fre*IU_free)./cx+BB*(Fre^2*IU_free)./((cx)^2);

% FF=FL;
% FS=-Sb_L*U_free_L(:,1);
% µ×²¿
% dxx=B_inf(:,2)-B_inf(1,2)*ones(size(B_inf,1),1);
% Bo_dx_mesh=diff(B_inf(:,2));
% Num_node=length(dxx);
% delay_EveryNode_t=dxx./cx;
% Txy=[Qxy(end,1);Qy(end,1)];
% U_bottom=U_free_L(end-1:end,1);
% for ii=1:Num_node
%     Sigmaxy(2*ii-1:2*ii,1)=Txy*exp(-1i*Fre*delay_EveryNode_t(ii));
%     B_u(2*ii-1:2*ii,1)    =U_bottom*exp(-1i*Fre*delay_EveryNode_t(ii));
% end
% NOCbottom=zeros(Num_node-1,2);
% NOCbottom(:,1)=1:Num_node-1; NOCbottom(:,2)=2:Num_node;
% Tao_xy=Sigmaxy(1:2:end,:);
% Sigma_y=Sigmaxy(2:2:end,:);
% Fx_bottom=zeros(Num_node,1);   Fy_bottom=zeros(Num_node,1);
% 
% for ele=1:Num_node-1
%     Fredom1=NOCbottom(ele, 1); Fredom2=NOCbottom(ele, 2);
%     Shape_bottom=[1/3*Bo_dx_mesh(ele),1/6*Bo_dx_mesh(ele);1/6*Bo_dx_mesh(ele),1/3*Bo_dx_mesh(ele)];
%     Fx_bottom(Fredom1:Fredom2,:)=Fx_bottom(Fredom1:Fredom2,:)+Shape_bottom*(Tao_xy(Fredom1:Fredom2,:));
%     Fy_bottom(Fredom1:Fredom2,:)=Fy_bottom(Fredom1:Fredom2,:)+Shape_bottom*(Sigma_y(Fredom1:Fredom2,:));
% end
% Force=zeros(2*Num_node,1);
% Force(1:2:end,1)=Fx_bottom;
% Force(2:2:end,1)=Fy_bottom;
% Feqb(degree_b,1)=+Force+stiffness_Half*B_u(:,1);
% Feq=Feql;%+Feqr;%+Feqb;
Uf=sparse(CoulpeK)\Feqr(:,1);
end

