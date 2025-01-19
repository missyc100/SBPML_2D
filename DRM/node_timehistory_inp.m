function ret = node_timehistory_inp(Free_Coorbb_value_x,Free_Coorbb_value_y,Free_Cooree_value_x,...
    Free_Cooree_value_y,Node_b_adjust,Node_e_adjust,KKbe,KKeb,dt,MMbe,MMeb,Free_Coorbb_value_ax,...
    Free_Coorbb_value_ay,Free_Cooree_value_ax,Free_Cooree_value_ay)
mkdir('output4')
% load load1.mat
% load KK.mat
% load node.mat
[m1,n1] = size(Free_Coorbb_value_x);
[m2,n2] = size(Free_Cooree_value_x);
ub = zeros(m1*2,n1);
ue = zeros(m2*2,n2);
ab = zeros(m1*2,n1);
ae = zeros(m2*2,n2);
for i=1:m1
    ub(i*2-1,:) = Free_Coorbb_value_x(i,:);
    ub(i*2,:) = Free_Coorbb_value_y(i,:);
    ab(i*2-1,:) = Free_Coorbb_value_ax(i,:);
    ab(i*2,:) = Free_Coorbb_value_ay(i,:);
end
for i=1:m2
    ue(i*2-1,:) = Free_Cooree_value_x(i,:);
    ue(i*2,:) = Free_Cooree_value_y(i,:);
    ae(i*2-1,:) = Free_Cooree_value_ax(i,:);
    ae(i*2,:) = Free_Cooree_value_ay(i,:);
end
Forceb=-MMbe*ae-KKbe*ue;
Forcee=MMeb*ab+KKeb*ub;
Force=[Forceb;Forcee];
parfor i109 = 1:m1
    ret = force_displace_gen(Forceb(i109*2-1,:), dt, Node_b_adjust(i109,1), 'x',0);
    ret = force_displace_gen(Forceb(i109*2,:), dt, Node_b_adjust(i109,1), 'y',0);
end
parfor i110 = 1:m2
    ret = force_displace_gen(Forcee(i110*2-1,:), dt, Node_e_adjust(i110,1), 'x',0);
    ret = force_displace_gen(Forcee(i110*2,:), dt, Node_e_adjust(i110,1), 'y',0);
end
ret = 'node_timehistory_inp finished'
% loading  = Force(1,:);
% node = 76;
% dt = 0.005;
% ret = loading_gen(loading, dt, node, 'x');

