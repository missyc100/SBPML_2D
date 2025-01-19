function [ABL,LABrs]=plus_Reig0(LB1,Bh,LABrs,BNL)
Sd=BNL;                          
ABh=zeros(Sd);

Degree_B=Sd;
ABh(Degree_B,Degree_B)=Bh;      
ABrs(Degree_B,Degree_B)=Bh/2;      
LABrs=LABrs+ABrs;
ABL=LB1+ABh;  
end

