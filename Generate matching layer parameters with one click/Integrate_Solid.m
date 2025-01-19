function [AArr,AArs,AAss,AMM] = Integrate_Solid(BNLR,NOC_BLR,BELR,B_infLR,nxLR)
%%%%%%≥ı ºæÿ’Û÷µ
BNL=BNLR;
AArr= zeros(BNL,BNL);  AArs= zeros(BNL,BNL); AAss= zeros(BNL,BNL);
AMM= zeros(BNL,BNL);
NOC_B=NOC_BLR;    BE=BELR;
for n = 1:BE
    [Arr,Ars,Ass,MM] = Ele_matrix(n,B_infLR,nxLR);
    degree = NOC_B(n,:);
    AArr(degree,degree) =  AArr(degree,degree) + Arr;
    AArs(degree,degree) =  AArs(degree,degree) + Ars;
    AAss(degree,degree) =  AAss(degree,degree) + Ass;
    AMM(degree,degree) =  AMM(degree,degree) + MM;
end
end