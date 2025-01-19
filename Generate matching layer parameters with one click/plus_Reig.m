function [Brr,Brs,Bss,AARS]=plus_Reig(Arr1,Ass_M,Ars_rsT,Ah,Ch,BNR,mmd)
Sd=BNR+mmd;

% AARS=zeros(Sd);
AAl=zeros(Sd);
ACl=zeros(Sd);

AAh=zeros(Sd);  
ACh=zeros(Sd);

% Degree_Up=1:mmu+1;
% AAu(Degree_Up,Degree_Up)=-Ahu;
% ACu(Degree_Up,Degree_Up)=-Chu;


Degree_layer=1:BNR;
AAl(Degree_layer,Degree_layer)=Arr1;
ACl(Degree_layer,Degree_layer)=Ass_M;
% AARS(Degree_layer,Degree_layer)=Ars1;

Degree_half=BNR:Sd;
AAh(Degree_half,Degree_half)=Ah;
ACh(Degree_half,Degree_half)=Ch;

Brs=zeros(Sd);
Brs(Degree_layer,Degree_layer)=1i*Ars_rsT;

% Ehh=Ars1+AEh;
Brr=AAl+AAh;
Bss=ACl+ACh;
end

