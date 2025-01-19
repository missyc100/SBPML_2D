z=1i*(5-3i)*3/0.5+4;
S=sqrt(1+z.^2);

JJ=200;
S2=2*(1+z.^2/2);
SS=zeros(size(z,2),JJ);
for jj=1:JJ
   S1=2*(1+z.^2/2)-z.^4./4./S2;
   S2=S1;
   SS(:,jj)= S1';
 end
S_C=1+z.^2/2-z.^4./4./S1;

plot(z,abs(S),'b-',z,abs(S_C),'r')
abs(SS);
