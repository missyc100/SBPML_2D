function [Fc,IU_free_B,Stress_free_Bxz,Stress_free_Byz,IU_free_R,Stress_free_Rxz,Stress_free_Ryz]...
    = incident_field(impedance,load_f,TimeB,Fre,TimeR,Free_Tx,Free_Ty,Free_ux,Ah,Ch,cx,mmd)
% ���䳡����
% ʱ���Թ���Fourier transformation 
% % ʵ�����ɶ�λ�����䳡λ��
IU_free_B=load_f.*exp(-1i*Fre*TimeB);
U_Free_B=Free_ux(end,1).*exp(-1i*Fre*TimeB);
IU_free_B(1:2,1)=U_Free_B(1:2,1);

Stress_free_Bxz=impedance(1,:)*1i*Fre*IU_free_B;
Stress_free_Byz=impedance(2,:)*1i*Fre*IU_free_B;

F_Free_Bx=Free_Tx(end,1).*exp(-1i*Fre*TimeB);
F_Free_By=Free_Ty(end,1).*exp(-1i*Fre*TimeB);

Stress_free_Bxz(1:2,1)=F_Free_Bx(1:2,1);
Stress_free_Byz(1:2,1)=F_Free_By(1:2,1);
% ���㸨�����ɶ�λ�����䳡λ��
IU_free_R=load_f.*exp(-1i*Fre*TimeR);
Stress_free_Rxz=impedance(1,:)*1i*Fre*IU_free_R;
Stress_free_Ryz=impedance(2,:)*1i*Fre*IU_free_R;

Co_M=Fre^2/cx^2*Ah+Ch;
De_R=1;De_A=2:mmd+1;


A11=Co_M(De_R,De_R);
A12=Co_M(De_R,De_A);
A21=Co_M(De_A,De_R);
A22=Co_M(De_A,De_A);
UA=-A22\A21*IU_free_R(end,1);

IU_free_R=[IU_free_R;UA];

Fc=-1i*Fre/cx*Ah*[IU_free_R(end,1);UA];
% IU_free_B=exp(-1i*Fre*TimeB)*Free_ux(end,1);
% Stress_free_Bxz=exp(-1i*Fre*TimeB)*Free_Tx(end,:);
% Stress_free_Byz=exp(-1i*Fre*TimeB)*Free_Ty(end,:);
% 
% 
% IU_free_R=Free_ux.*exp(-1i*Fre*TimeR);
% Stress_free_Rxz=Free_Tx.*exp(-1i*Fre*TimeR);
% Stress_free_Ryz=Free_Ty.*exp(-1i*Fre*TimeR);
% ��⸨�����ɶȴ�λ��




end

