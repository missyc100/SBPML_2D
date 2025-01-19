function [Ahs,Bhs,Chs] = Augmented_matrix(CsH,GH)

Ahs=GH*CsH/2;
Bhs=GH/1i;
Chs=GH/CsH;
end

