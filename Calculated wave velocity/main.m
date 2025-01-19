close all;clear all;
clc
rs = 2000;% density
v = 0.3;% Poisson's ratio
Cs = 200;% Shear wave velocity
G = Cs^2*rs;% Shear modulus
E = G*(2*(1+v));% Modulus of elasticity
lame = (2*G*v)./(1-2*v);% lame constant
Cp = sqrt((lame+2*G)./rs);% Compression wave velocity
%% elements size
f = 9;% Incident wave frequency
L_max = Cp/f;% Maximum wavelength
L = Cs/f;% Minimum wavelength
x = 1/10*L;% 
%% the length of time step
t = x/Cp;
t2 = 1/15/f;

m = 2;% The order of the stretching function
b_ele = 2*10;% 10 times the recommended elements size
PHY_LPML = 20;% SBPML domain thickness
R = 10^-9;% Absorption coefficient
C_ref = Cs;% Shear wave velocity

alfa0=(m+1)*b_ele/(2*(PHY_LPML))*log(1/abs(R));% Attenuation coefficient of waves
beta0=(m+1)*C_ref/(2*(PHY_LPML))*log(1/abs(R));% Space tensile coefficient



