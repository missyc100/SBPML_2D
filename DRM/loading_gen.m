function ret = loading_gen(loading, dt, node, dime)
% load load1.mat
% loading  = Free_Coorbb_value_x(1,:);
% node = 1;
title1 = '*Amplitude';
title2 = ' name=Amp-';
% title = char(title);
space1 = '             ';
% dt=0.005;
ntime = length(loading);
t=(0:dt:(ntime-1)*dt);
% output_x = zeros(3,3);
output_x_initial{1,1} = title1;
output_x_initial{1,2} = [title2,char(string(node)),dime];
writecell(output_x_initial,['output2/output_',char(string(node)),dime,'.txt'])
for i = 1:ntime
    mod1 = mod(i,4);
    row1 = fix(i/4);
    if mod1 == 0
        mod1 = 4;
        output_x{row1,mod1*2-1} = char(string(t(i)));
        output_x{row1,mod1*2} = char(string(loading(i)));
    else
        output_x{row1+1,mod1*2-1} = char(string(t(i)));
        output_x{row1+1,mod1*2} = char(string(loading(i)));
    end
end
output_x_middel = output_x(1:end-1,:);
writecell(output_x_middel,['output2/output_',char(string(node)),dime,'.txt'],'WriteMode','append')
for i=1:mod1
    output_x_end{1,i*2-1} = output_x{row1+1,i*2-1};
    output_x_end{1,i*2} = output_x{row1+1,i*2};
end
writecell(output_x_end,['output2/output_',char(string(node)),dime,'.txt'],'WriteMode','append')
ret = 'finish';




% writecell(rgb,'C.xls','WriteMode','append')
% A = readcell('output_x.txt');
% A = readcell('test1.txt');
% A{1,3} = ' ';
% A{1,4} = ' ';
% A{1,5} = ' ';
% A{1,6} = ' ';
% A{1,7} = ' ';
% A{1,8} = ' ';
% writecell(A,'output_x.txt')