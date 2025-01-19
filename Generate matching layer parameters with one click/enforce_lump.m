function [de,pressure] = enforce_lump(loading,XofN_R,Force)
XofNA=XofN_R;
% esp=1.0E-8;
% Force= importdata('Force.txt');
% Force= Force(find( ~ isnan(Force)));

if size(Force,1)==1
   Force=Force.';
end
Force_N=[Force,XofN_R(Force,:)];
Force_N=sortrows(Force_N, [2,3]);

% ro_imposeY=find(abs(XofNA(:,2)==0 );  
% ro_imposeX=find(abs(XofNA(ro_imposeY,1)-0)<=0.4+0.01);  
% X_co=XofNA(ro_imposeY(ro_imposeX),:);
% B_surface=[ro_imposeY(ro_imposeX),X_co];
% B_surface1=sortrows(B_surface,2);

CoorXY=Force_N(:,1);
pressure=loading;
de=zeros(size(CoorXY,1),1);
% de(1:2:end,:)=2*CoorXY-1;
de(1:end,:)=2*CoorXY-1;

% 固定边界条件
% X_max=max(XofNA(:,1));X_min=min(XofNA(:,1));
% Y_max=max(XofNA(:,2));
% % P=[X_min,X_max,Y_max].';
% [ro1,co]=find(abs(XofNA(:,1)-X_min)<esp);
% [ro2,co]=find(abs(XofNA(:,1)-X_max)<esp);
% [ro3,co]=find(abs(XofNA(:,2)-Y_max)<esp);
% N_B=unique([ro1;ro2;ro3]);

% N_B= importdata('Fixed_boundary.txt');
% N_B= N_B(find( ~ isnan(N_B)));
end

