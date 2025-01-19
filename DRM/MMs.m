   function [ME] = MMs(n,XofNs,NOCs,cx,lamds,Gs,rs)%%¸úFreÃ»¹ØÏµ
Y=XofNs(:,2);
NOC_B=NOCs(:,2:3);
ME=zeros (4,4);
Y1= Y(NOC_B(n, 1),1);
Y2= Y(NOC_B(n, 2),1);
dy=-(Y2-Y1); 
ME(1,1)=-2*(lamds(n)+2*Gs(n))/cx^2+2*rs(n);   ME(1,2)=0;                        ME(1,3)=-(lamds(n)+2*Gs(n))/cx^2+rs(n);       ME(1,4)=0;                    
ME(2,1)=0;                                   ME(2,2)=-2*Gs(n)/cx^2+2*rs(n);    ME(2,3)=0;                                   ME(2,4)=-Gs(n)/cx^2+rs(n);                    
ME(3,1)=-(lamds(n)+2*Gs(n))/cx^2+rs(n);       ME(3,2)=0;                        ME(3,3)=-2*(lamds(n)+2*Gs(n))/cx^2+2*rs(n);   ME(3,4)=0;                    
ME(4,1)=0;                                   ME(4,2)=-Gs(n)/cx^2+rs(n);        ME(4,3)=0;                                   ME(4,4)=-2*Gs(n)/cx^2+2*rs(n);
ME=ME*dy/6;   
end

