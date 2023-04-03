%% Dynamics of Safety Embedded Planar Double Integrator
function [f_z,fx_z,fu_z,z0,zf]=BaS_dyn_symb(o11,o12,r1,o21,o22,r2,o31,o32,r3,o41,o42,r4,gamma,x0,xf);
syms x1 x2 x3 x4 u1 u2
x=[x1;x2;x3;x4]; u=[u1;u2];

% define safety paramters:
h1= (x1-o11)^2+(x2-o12)^2-r1^2; h1x=[2*(x1-o11),2*(x2-o12),0,0];
h2= (x1-o21)^2+(x2-o22)^2-r2^2; h2x=[2*(x1-o21),2*(x2-o22),0,0];
h3= (x1-o31)^2+(x2-o32)^2-r3^2; h3x=[2*(x1-o31),2*(x2-o32),0,0];
h4= (x1-o41)^2+(x2-o42)^2-r4^2; h4x=[2*(x1-o41),2*(x2-o42),0,0];
%define barrier state
h1_orig=o11^2+o12^2-r1^2; c1=1/h1_orig; 
h2_orig=o21^2+o22^2-r2^2; c2=1/h2_orig; 
h3_orig=o31^2+o32^2-r3^2; c3=1/h3_orig; 
h4_orig=o41^2+o42^2-r4^2; c4=1/h4_orig; 
c=c1+c2+c3+c4;

syms x5;
zc=x5+c;
% call system's dynamics (maybe symbolic)
f=PlanDInt_dynamics(x,u,false);

h=[h1;h2;h3;h4]; hx=[h1x;h2x;h3x;h4x]; f_z=0;
for ii=1:length(h);
    sum_Bofh=0;
    for jj=1:length(h)
        if jj~=ii
            Bofh=1/(h(jj));
        sum_Bofh=sum_Bofh+Bofh;
        end
    end
    phi1=-(zc-sum_Bofh)^2*h(ii)-(zc-sum_Bofh);
    phi0=-(zc-sum_Bofh)^2;
f_z= f_z +(-gamma*phi1 + phi0 * hx(ii,:)*f);
end


% f_z=- gamma*(h*zc^2-zc)- zc^2*(hx*f);
% hkplus=(x3-o1)^2 + (x4-o2)^2-r^2; 
% f_z= 1/hkplus - c;
fx_z=jacobian(f_z,[x;x5]);
fu_z=jacobian(f_z,u); 


h10=double(subs(h1,[x1,x2,x3,x4],x0')); h30=double(subs(h3,[x1,x2,x3,x4],x0')); 
h20=double(subs(h2,[x1,x2,x3,x4],x0')); h40=double(subs(h4,[x1,x2,x3,x4],x0')); 
z0=1/h10+1/h20+1/h30+1/h40-c; 
h1f=double(subs(h1,[x1,x2,x3,x4],xf')); h3f=double(subs(h3,[x1,x2,x3,x4],xf')); 
h2f=double(subs(h2,[x1,x2,x3,x4],xf')); h4f=double(subs(h4,[x1,x2,x3,x4],xf')); 
zf=1/h1f+1/h2f+1/h3f+1/h4f-c; 
end
