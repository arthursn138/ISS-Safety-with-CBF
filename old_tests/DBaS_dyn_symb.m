%% Dynamics of Safety Embedded Planar Double Integrator
function [f_z,fx_z,fu_z,z0,zf]=DBaS_dyn_symb(x0,xf,f_dyn,safe_par);
o11=safe_par.o11; o12=safe_par.o12; r1=safe_par.r1;
o21=safe_par.o21; o22=safe_par.o22; r2=safe_par.r2;
o31=safe_par.o31; o32=safe_par.o32; r3=safe_par.r3;
o41=safe_par.o41; o42=safe_par.o42; r4=safe_par.r4;
syms x1 x2 x3 x4 u1 u2
x=[x1;x2;x3;x4]; u=[u1;u2];
xplus=f_dyn(x,u); x1plus=xplus(1); x2plus=xplus(2);
% define safety equations:
h1= @(x1,x2) (x1-o11)^2+(x2-o12)^2-r1^2; 
h2= @(x1,x2) (x1-o21)^2+(x2-o22)^2-r2^2; 
h3= @(x1,x2) (x1-o31)^2+(x2-o32)^2-r3^2; 
h4= @(x1,x2) (x1-o41)^2+(x2-o42)^2-r4^2; 
h1plus= (x1plus-o11)^2+(x2plus-o12)^2-r1^2;
h2plus= (x1plus-o21)^2+(x2plus-o22)^2-r2^2;
h3plus= (x1plus-o31)^2+(x2plus-o32)^2-r3^2;
h4plus= (x1plus-o41)^2+(x2plus-o42)^2-r4^2;
%define barrier state
% beta = 1/h1 + 1/h2 + 1/h3 + 1/h4 + 1/h5;
beta_plus= 1/h1plus + 1/h2plus + + 1/h3plus + 1/h4plus;
% beta_x= b1x + b2x + b3x + b4x+ b5x;
% beta_xx=b1xx + b2xx + b3xx + b4xx + b5xx;
% det:
f_z= beta_plus; %det barrier
syms x5
fx_z=jacobian(f_z,[x;x5]);
fu_z=jacobian(f_z,u); 

z0=1/h1(x0(1),x0(2)) +1/h2(x0(1),x0(2)) +1/h3(x0(1),x0(2)) +1/h4(x0(1),x0(2));
zf=1/h1(xf(1),xf(2)) +1/h2(xf(1),xf(2)) +1/h3(xf(1),xf(2)) +1/h4(xf(1),xf(2));

% symbolic to function handle
f_z=matlabFunction(f_z);
fx_z=matlabFunction(fx_z);
fu_z=matlabFunction(fu_z);

% h.h1=matlabFunction(h1);
% h.h2=matlabFunction(h2);
% h.h3=matlabFunction(h3);
% h.h4=matlabFunction(h4);
% h.h5=matlabFunction(h5);
end
