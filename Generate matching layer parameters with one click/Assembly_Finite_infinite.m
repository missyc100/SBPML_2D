function [ MM,CC,KK ]=Assembly_Finite_infinite(Fi_MM,Fi_KK,MM1,CC1,KK1,SN,N_reg)
MM_r=zeros(SN);  MM_pml=MM_r; 
CC_pml=MM_r;  
KK_r=MM_r;       KK_pml=MM_r;


De_r=1:size(Fi_MM,1);


Fi_N=2*size(N_reg,1);
De_pml=Fi_N+1:SN;

MM_r(De_r,De_r)=Fi_MM;
MM_pml(De_pml,De_pml)=MM1;

CC_pml(De_pml,De_pml)=CC1;

KK_r(De_r,De_r)=Fi_KK;
KK_pml(De_pml,De_pml)=KK1;

MM=MM_r+MM_pml;
CC=CC_pml;
KK=KK_r+KK_pml;
clear MM_r MM_pml CC_pml KK_r KK_pml

end


