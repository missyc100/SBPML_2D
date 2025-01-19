function [A1,B1,C1,D1] = Dry_Co( ABrr_R,ABrs_R,ABss_R,AMM_R,Fre)
A1 = ABrr_R;
B1 = ABrs_R-ABrs_R';
C1 = ABss_R-Fre^2*AMM_R;
D1 = ABrs_R;
end

