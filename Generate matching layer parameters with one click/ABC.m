function [S_abc] = ABC(lamdN, GsN,rs,mp,np,ms,ns,cx,PorSV)
%UNTITLED 此处显示有关此函数的摘要
D_elastic=[lamdN+2*GsN,lamdN,0;lamdN,lamdN+2*GsN,0;0,0 GsN];
Strain=[mp^2,-ns*ms; np^2, ns*ms;2*mp*np, ms^2-ns^2 ];
D_d=-cx*[mp^2,-ms*ns;np*mp,ms^2];
S_ABC=D_elastic*Strain/D_d;
S_abc=-[S_ABC(3,:);S_ABC(2,:)];

if cx==Inf && PorSV==1
    S_abc=[0,0;0,sqrt((lamdN+2*GsN)*rs)];
end
if cx==Inf && PorSV==2
    S_abc=[sqrt((GsN)*rs),0,;0,0];
end
end

