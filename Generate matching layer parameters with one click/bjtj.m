function [UU]=bjtj(U_free_L,L_B1_inf,L_B2_inf,nL1,nL2,...
    U_free_R,R_B1_inf,R_B2_inf,nR1,nR2,...
    Sigma1_LSx,Sigma1_LSy,Sigma1_LSxy, Sigma1_RSx,Sigma1_RSy,Sigma1_RSxy,...
    Sigma2_LSx,Sigma2_LSy,Sigma2_LSxy, Sigma2_RSx,Sigma2_RSy,Sigma2_RSxy,...
    Sigma2_LDpx,Sigma2_LDpy,Sigma2_RDpx,Sigma2_RDpy,...
    CoulpeK,Sb_L_A,Sb_R_A,plus_BLL,plus_BRR,Bouconditon)
Num=size(CoulpeK,1);
Feql=zeros(Num,1); Feqr=zeros(Num,1); %Feqb=zeros(2*NNs,1);
%% Dry Soil
% Left
L_Num_node=size(L_B1_inf,1);
NOCL=zeros(L_Num_node-1,2);
NOCL(:,1)=1:L_Num_node-1; NOCL(:,2)=2:L_Num_node;
L_dy_mesh0=diff(L_B1_inf(:,2:3));
L_dy_mesh=sqrt(L_dy_mesh0(:,1).^2+L_dy_mesh0(:,2).^2);
Fx_L1=zeros(L_Num_node,1); Fy_L1=zeros(L_Num_node,1);
for ele=1:L_Num_node-1
    Fredom1=NOCL(ele, 1); Fredom2=NOCL(ele, 2);
    Shape_L=[1/3*L_dy_mesh(ele),1/6*L_dy_mesh(ele);1/6*L_dy_mesh(ele),1/3*L_dy_mesh(ele)];
    Fx_L1(Fredom1:Fredom2,:)=Fx_L1(Fredom1:Fredom2,:)+Shape_L*(-nL1(ele,2)*Sigma1_LSx(Fredom1:Fredom2,:)+nL1(ele,1)*Sigma1_LSxy(Fredom1:Fredom2,:));
    Fy_L1(Fredom1:Fredom2,:)=Fy_L1(Fredom1:Fredom2,:)+Shape_L*(nL1(ele,1)*Sigma1_LSy(Fredom1:Fredom2,:)-nL1(ele,2)*Sigma1_LSxy(Fredom1:Fredom2,:));
end
FL1=zeros(2*L_Num_node,1);
FL1(1:2:end,:)=Fx_L1;FL1(2:2:end,:)=Fy_L1;
% Right
R_Num_node=size(R_B1_inf,1);
NOCR=zeros(R_Num_node-1,2);
NOCR(:,1)=1:R_Num_node-1; NOCR(:,2)=2:R_Num_node;
R_dy_mesh0=diff(R_B1_inf(:,2:3));
R_dy_mesh=sqrt(R_dy_mesh0(:,1).^2+R_dy_mesh0(:,2).^2);
Fx_R1=zeros(R_Num_node,1); Fy_R1=zeros(R_Num_node,1);
for ele=1:R_Num_node-1
    Fredom1=NOCR(ele, 1); Fredom2=NOCR(ele, 2);
    Shape_R=[1/3*R_dy_mesh(ele),1/6*R_dy_mesh(ele);1/6*R_dy_mesh(ele),1/3*R_dy_mesh(ele)];
    Fx_R1(Fredom1:Fredom2,:)=Fx_R1(Fredom1:Fredom2,:)+Shape_R*(nR1(ele,2)*Sigma1_RSx(Fredom1:Fredom2,:)-nR1(ele,1)*Sigma1_RSxy(Fredom1:Fredom2,:));
    Fy_R1(Fredom1:Fredom2,:)=Fy_R1(Fredom1:Fredom2,:)+Shape_R*(-nR1(ele,1)*Sigma1_RSy(Fredom1:Fredom2,:)+nR1(ele,2)*Sigma1_RSxy(Fredom1:Fredom2,:));
end
FR1=zeros(2*R_Num_node,1);
FR1(1:2:end)=Fx_R1;FR1(2:2:end,:)=Fy_R1;
%% Saturated Soil
% Left
L_Num_node2=size(L_B2_inf,1);
NOCL2=zeros(L_Num_node2-1,2);
NOCL2(:,1)=1:L_Num_node2-1; NOCL2(:,2)=2:L_Num_node2;
L_dy_mesh02=diff(L_B2_inf(:,2:3));
L_dy_mesh2=sqrt(L_dy_mesh02(:,1).^2+L_dy_mesh02(:,2).^2);
Fx_L2=zeros(L_Num_node2,1); Fy_L2=zeros(L_Num_node2,1); Fp_L2=zeros(L_Num_node2,1);
for ele=1:L_Num_node2-1
    Fredom1=NOCL2(ele, 1); Fredom2=NOCL2(ele, 2);
    Shape_L2=[1/3*L_dy_mesh2(ele),1/6*L_dy_mesh2(ele);1/6*L_dy_mesh2(ele),1/3*L_dy_mesh2(ele)];
    Fx_L2(Fredom1:Fredom2,:)=Fx_L2(Fredom1:Fredom2,:)+Shape_L2*(-nL2(ele,2)*Sigma2_LSx(Fredom1:Fredom2,:)+nL2(ele,1)*Sigma2_LSxy(Fredom1:Fredom2,:));
    Fy_L2(Fredom1:Fredom2,:)=Fy_L2(Fredom1:Fredom2,:)+Shape_L2*(nL2(ele,1)*Sigma2_LSy(Fredom1:Fredom2,:)-nL2(ele,2)*Sigma2_LSxy(Fredom1:Fredom2,:));
    Fp_L2(Fredom1:Fredom2,:)=Fp_L2(Fredom1:Fredom2,:)+Shape_L2*(-nL2(ele,2)*Sigma2_LDpx(Fredom1:Fredom2,:)+nL2(ele,1)*Sigma2_LDpy(Fredom1:Fredom2,:));
end
FL2=zeros(2*L_Num_node2,1);
FL2(1:2:end)=Fx_L2;FL2(2:2:end,:)=Fy_L2;
% Right
R_Num_node2=size(R_B2_inf,1);
NOCR2=zeros(R_Num_node2-1,2);
NOCR2(:,1)=1:R_Num_node2-1; NOCR2(:,2)=2:R_Num_node2;
R_dy_mesh02=diff(R_B2_inf(:,2:3));
R_dy_mesh2=sqrt(R_dy_mesh02(:,1).^2+R_dy_mesh02(:,2).^2);
Fx_R2=zeros(R_Num_node2,1); Fy_R2=zeros(R_Num_node2,1); Fp_R2=zeros(R_Num_node2,1);
for ele=1:R_Num_node2-1
    Fredom1=NOCR2(ele, 1); Fredom2=NOCR2(ele, 2);
    Shape_R2=[1/3*R_dy_mesh2(ele),1/6*R_dy_mesh2(ele);1/6*R_dy_mesh2(ele),1/3*R_dy_mesh2(ele)];
    Fx_R2(Fredom1:Fredom2,:)=Fx_R2(Fredom1:Fredom2,:)+Shape_R2*(nR2(ele,2)*Sigma2_RSx(Fredom1:Fredom2,:)-nR2(ele,1)*Sigma2_RSxy(Fredom1:Fredom2,:));
    Fy_R2(Fredom1:Fredom2,:)=Fy_R2(Fredom1:Fredom2,:)+Shape_R2*(-nR2(ele,1)*Sigma2_RSy(Fredom1:Fredom2,:)+nR2(ele,2)*Sigma2_RSxy(Fredom1:Fredom2,:));
    Fp_R2(Fredom1:Fredom2,:)=Fp_R2(Fredom1:Fredom2,:)+Shape_R2*(nR2(ele,2)*Sigma2_RDpx(Fredom1:Fredom2,:)-nR2(ele,1)*Sigma2_RDpy(Fredom1:Fredom2,:));
end
FR2=zeros(2*R_Num_node2,1);
FR2(1:2:end)=Fx_R2;FR2(2:2:end,:)=Fy_R2;

% Commom nodal
L_Num=2*(L_Num_node+L_Num_node2-1)+L_Num_node2;
FLD=zeros(L_Num,1); FLS=zeros(L_Num,1); FLP=zeros(L_Num,1);
FLD(1:2*L_Num_node,1)=FL1;
FLS(2*L_Num_node-1:2*(L_Num_node+L_Num_node2-1),1)=FL2;
FLP(2*(L_Num_node+L_Num_node2-1)+1:L_Num,:)=Fp_L2;
FL=FLD+FLS+FLP;
% Right
R_Num=2*(R_Num_node+R_Num_node2-1)+R_Num_node2;
FRD=zeros(R_Num,1); FRS=zeros(R_Num,1); FRP=zeros(R_Num,1);
FRD(1:2*R_Num_node,1)=FR1;
FRS(2*R_Num_node-1:2*(R_Num_node+R_Num_node2-1),1)=FR2;
FRP(2*(R_Num_node+R_Num_node2-1)+1:R_Num,1)=Fp_R2;
FR=FRD+FRS+FRP;

Feqr(plus_BRR,1)=FR(:,1)+Sb_R_A*U_free_R(:,1);
Feql(plus_BLL,1)=FL(:,1)-Sb_L_A*U_free_L(:,1);
% µ×²¿ans
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
Feq=Feql+Feqr;%+Feqb;
CoulpeK(Bouconditon,:)=0;
CoulpeK(:,Bouconditon)=0;
CoulpeK(Bouconditon,Bouconditon)=eye(length(Bouconditon));
Feq(Bouconditon,1)=0;
UU=sparse(CoulpeK)\Feq(:,1);
end

