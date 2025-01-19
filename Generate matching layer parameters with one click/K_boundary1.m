function  [CoulpeK,Plus_BRR]=K_boundary1(MM,KK,Sb_R,plus_BR,Fre,mmd,plus_BB,Absorb_matrices,plus_BL,Sb_L)
Keq=-Fre^2*MM+KK;
Real_N=size(Keq,1);
Freedom=Real_N+mmd;
Plus_BRR=[plus_BR;[Real_N+1:Freedom]'];
% Plus_BRR=[plus_BR];

K_boundaryR=zeros(Freedom);
K_boundaryL=zeros(Freedom);
K_boundaryB=zeros(Freedom);

Fi_K=zeros(Freedom);
% % K_boundaryL=zeros(size(Keq,1));
K_boundaryL(plus_BL,plus_BL)=Sb_L; 
K_boundaryR(Plus_BRR,Plus_BRR)=Sb_R;
K_boundaryB(plus_BB(1:end,1),plus_BB(1:end,1))=1i*Fre*Absorb_matrices(1:end,1:end);
Fi_K(1:Real_N,1:Real_N)=Keq;
CoulpeK=K_boundaryR+Fi_K+K_boundaryB+K_boundaryL;

end

