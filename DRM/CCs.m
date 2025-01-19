function [Cs] = CCs(n,XofNS,NOCS,cx,lamds,Gs,rs)
Cs=zeros (4,4);
Cs(1,1)=0;                     Cs(1,2)=-(lamds(n)-Gs(n))/2/cx;     Cs(1,3)=0;                    Cs(1,4)=(lamds(n)+Gs(n))/2/cx;       
Cs(2,1)=(lamds(n)-Gs(n))/2/cx;  Cs(2,2)=0;                         Cs(2,3)=(lamds(n)+Gs(n))/2/cx; Cs(2,4)=0;        
Cs(3,1)=0;                     Cs(3,2)=(-lamds(n)-Gs(n))/2/cx;     Cs(3,3)=0;                    Cs(3,4)=(lamds(n)-Gs(n))/2/cx;        
Cs(4,1)=(-lamds(n)-Gs(n))/2/cx; Cs(4,2)=0;                         Cs(4,3)=(-lamds(n)+Gs(n))/2/cx;Cs(4,4)=0;        
end 

                                                                                                                                                                          