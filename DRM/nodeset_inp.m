function ret = nodeset_inp(Node_b_adjust,Node_e_adjust,ampx,ampy)
mkdir('output1')
% load node.mat
Node = [Node_b_adjust;Node_e_adjust];
m = length(Node);
% *Cload, amplitude=Amp-1
% Set-4, 2, -1.
for i=1:m
    output1{1,1} = '*Nset';
    output1{1,2} = [' nset=Set-',char(string(i+10))];
    output1{1,3} = ' instance=Part-1-1';
    writecell(output1,'output1/output.txt','WriteMode','append')
    output2{1,1} = Node(i);
    writecell(output2,'output1/output.txt','WriteMode','append')

    output3_0{1,1} = ['** Name: Load-',char(string(i*2-1)),'   Type: Concentrated force'];
    writecell(output3_0,'output1/output2.txt','WriteMode','append')
    output3{1,1} = '*Cload';
    output3{1,2} = ['amplitude=Amp-',char(string(Node(i))),'x'];
    writecell(output3,'output1/output2.txt','WriteMode','append')
    output4{1,1} = ['Set-',char(string(i+10))];
    output4{1,2} = 1;
    output4{1,3} = ampx;
    writecell(output4,'output1/output2.txt','WriteMode','append')

    output5_0{1,1} = ['** Name: Load-',char(string(i*2)),'   Type: Concentrated force'];
    writecell(output5_0,'output1/output2.txt','WriteMode','append')
    output5{1,1} = '*Cload';
    output5{1,2} = ['amplitude=Amp-',char(string(Node(i))),'y'];
    writecell(output5,'output1/output2.txt','WriteMode','append')
    output6{1,1} = ['Set-',char(string(i+10))];
    output6{1,2} = 2;
    output6{1,3} = ampy;
    writecell(output6,'output1/output2.txt','WriteMode','append')
    %     output4{1,1} = '*Cload';
    %     output4{1,2} = ['*amplitude=Amp-',Node(i),'y'];
    
end
ret = 'nodeset_inp finished'

