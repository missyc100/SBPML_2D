function ret = force_displace_gen(displace, dt, node, dime, nndt)
title1 = '*Amplitude';
title2 = ' name=Amp-';
space1 = '             ';
ntime = length(displace);
t=(0:dt:(ntime-1)*dt);
if nndt ~= 0
    t = t+nndt;
    t = [0,t];
    displace = [0,displace];
    ntime = ntime+1;
end
output_x_initial{1,1} = title1;
output_x_initial{1,2} = [title2,char(string(node)),dime];
dir = ['output4/outputP',char(string(node)),char(dime),'.txt'];
% dir = ['output4/output-P','.txt'];
writecell(output_x_initial,dir,'WriteMode','append')
for i = 1:ntime
    mod1 = mod(i,4);
    row1 = fix(i/4);
    if mod1 == 0
        mod1 = 4;
        output_x{row1,mod1*2-1} = char(num2str(t(i),6));
        output_x{row1,mod1*2} = char(num2str(displace(i),'%.20f'));
    else
        output_x{row1+1,mod1*2-1} = char(num2str(t(i),15));
        output_x{row1+1,mod1*2} = char(num2str(displace(i),'%.20f'));
    end
end
output_x_middel = output_x(1:end-1,:);
writecell(output_x_middel,dir,'WriteMode','append')
if mod1~=4
    for i=1:mod1
        output_x_end{1,i*2-1} = output_x{row1+1,i*2-1};
        output_x_end{1,i*2} = output_x{row1+1,i*2};
    end
    writecell(output_x_end,dir,'WriteMode','append')
end

ret = 'displace_gen_com_finish';