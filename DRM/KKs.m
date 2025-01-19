function [Ks] = KKs(n,XofNs,NOCs,cx,lamds,Gs,rs)
NOC_B=NOCs(:,2:3);
Y=XofNs(:,2);
Y1= Y(NOC_B(n, 1),1);
Y2= Y(NOC_B(n, 2),1);
dy=-(Y2-Y1);
Ks=zeros (4,4);
Ks(1,1)=Gs(n)/dy; Ks(1,2)=0;                    Ks(1,3)=-Gs(n)/dy; Ks(1,4)=0                 ;
Ks(2,1)=0;        Ks(2,2)=(lamds(n)+2*Gs(n))/dy; Ks(2,3)=0;         Ks(2,4)=-(lamds(n)+2*Gs(n))/dy ;
Ks(3,1)=-Gs(n)/dy;Ks(3,2)=0;                    Ks(3,3)=Gs(n)/dy;  Ks(3,4)=0                 ;
Ks(4,1)=0;        Ks(4,2)=-(lamds(n)+2*Gs(n))/dy;Ks(4,3)=0;         Ks(4,4)=(lamds(n)+2*Gs(n))/dy;
end 

                                                                                                                                                                          