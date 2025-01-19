function [ BB,ARR, ARS ] = Stiffness_Co(Ars1,Arr1,Brr_half,Ars_half,Arr_half,BNR,mmd )

Sd=BNR+mmd;

ARS1=zeros(Sd);
ARS2=zeros(Sd);
ARS3=zeros(Sd);


ARR1=zeros(Sd);
ARR2=zeros(Sd);
ARR3=zeros(Sd);
BB2=zeros(Sd);


% Degree_Up=1:mmu+1;
% ARS3(Degree_Up,Degree_Up)=Ars_up;
% ARR3(Degree_Up,Degree_Up)=Arr_up;
% BB3(Degree_Up,Degree_Up)=Brr_up;




Degree_layer=1:BNR;
ARS1(Degree_layer,Degree_layer)=Ars1;
ARR1(Degree_layer,Degree_layer)=Arr1;

Degree_half=Sd-mmd:Sd;
ARS2(Degree_half,Degree_half)=Ars_half;
ARR2(Degree_half,Degree_half)=Arr_half;
BB2(Degree_half,Degree_half)=Brr_half;



ARS=ARS3+ARS2+ARS1;
ARR=ARR3+ARR2+ARR1;
BB=BB2;




end

