function [MM,CC,KK] = integrate_S_1D(BNLR,...
    NOC_BLR,BELR,B_infLR,Ah,cx)
BNL=BNLR;
AA11= zeros(BNL,BNL); AA22= zeros(BNL,BNL);
AMM= zeros(BNL,BNL);
NOC_B=NOC_BLR;    BE=BELR;
for n = 1:BE
    [A11,A22,MM] = Ele_matrix_1D(n,B_infLR);
    degree = NOC_B(n,:);
    AA11(degree,degree) =  AA11(degree,degree) + A11;
    AA22(degree,degree) =  AA22(degree,degree) + A22;
    AMM(degree,degree) =  AMM(degree,degree) + MM;
end
  MM=(AMM-AA11./cx^2);
  CC=zeros(BNL,BNL);
  CC(BNL,BNL)=Ah;
  KK=AA22;
end