function [KKbe,KKeb,MMbe,MMeb]=GENERATE_KK(Node_b_adjust,Node_e_adjust,FileName1,FileName2)
%% load matlab.mat
for i = 1:length(Node_b_adjust)
    dof_b(i*2-1,1) = Node_b_adjust(i)*2-1;
    dof_b(i*2,1) = Node_b_adjust(i)*2;
end
for i = 1:length(Node_e_adjust)
    dof_e(i*2-1,1) = Node_e_adjust(i)*2-1;
    dof_e(i*2,1) = Node_e_adjust(i)*2;
end
% C = mtx2Str(FileName1);

C = mtx2Str(FileName1);
D = mtx2Double(C);
if size(D,2) == 6
    D(:,6)=[];
end
save MTX\KK.mtx -ascii D


C = [];
D = [];
C = mtx2Str(FileName2);
D = mtx2Double(C);
if size(D,2) == 6
    D(:,6)=[];
end
save MTX\MM.mtx -ascii D
%% KK
dofsPerNode = 2;
[Mout1,out1] = readMtx('MTX\KK.mtx','matrix',...
'matrix input','global',dofsPerNode);
% [Mout1,out1] = readMtx('Job-12-KKMM_STIF2.mtx','matrix',...
% 'matrix input','global',dofsPerNode);
% load MTX/KK.mat
% K = KK;
K = Mout1{1,1};
% TK = K';
% TK(logical(speye(length(K))))=0;
% M = full(M);
% for i = 1:length(M)
%     for j = 1:length(M)
%         if i==j
%             TM(i,j) = 0;
%         end
%     end
% end
% K = K+TK;
KKbe = K(dof_b,dof_e);
KKbe = full(KKbe);
KKeb = K(dof_e,dof_b);
KKeb = full(KKeb);
% KKbb = K(dof_b,dof_b);
% KKee = K(dof_e,dof_e);
%% MM
[Mout2,out2] = readMtx('MTX\MM.mtx','matrix',...
'matrix input','global',dofsPerNode);
% [Mout1,out1] = readMtx('Job-12-1_STIF2.mtx','matrix',...
% 'matrix input','global',dofsPerNode);
M = Mout2{1,1};
TM = M';
TM(logical(speye(length(M))))=0;
% M = full(M);
% for i = 1:length(M)
%     for j = 1:length(M)
%         if i==j
%             TM(i,j) = 0;
%         end
%     end
% end
M = M+TM;
MMbe = M(dof_b,dof_e);
MMbe = full(MMbe);
MMeb = M(dof_e,dof_b);
MMeb = full(MMeb);
MMbb = M(dof_b,dof_b);
MMee = M(dof_e,dof_e);
%%
ret = 'GENERATE_KK_MM finished'
%%






