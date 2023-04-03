%% Dynamics of Safety Embedded Planar Double Integrator
function [fbar,fbarx,fbaru]=SafeEmbd_PlanDInt_dynamics(x,u,deriv_bool,f_z,fx_z,fu_z,dt)
% syms x1 x2 x3 x4 x5 u1 u2

x1=x(1);x2=x(2);x3=x(3);x4=x(4);x5=x(5);
% define safety paramters:
% f_z=double(subs(f_z,[x1,x2,x3,x4,x5,u1,u2],[x',u']));
% f_z=x5/10 - (x3*(2*x1 - 4) + x4*(2*x2 - 4))*(x5 + 16/127)^2 - ((x5 + 16/127)^2*((x1 - 2)^2 + (x2 - 2)^2 - 1/16))/10 + 8/635;

f_z=f_z(x1,x2,x3,x4);
f=[x(1);x(2);x(3);x(4)]+dt*[x(3); x(4); u(1); u(2)];
fbar=[f;f_z];
    if deriv_bool
        % dynamics gradients
%         fx_z=double(subs(fx_z,[x1,x2,x3,x4,x5,u1,u2],[x',u']));
%         fx_z=[ - 2*x3*(x5 + 16/127)^2 - ((2*x1 - 4)*(x5 + 16/127)^2)/10, - 2*x4*(x5 + 16/127)^2 - ((2*x2 - 4)*(x5 + 16/127)^2)/10, -(2*x1 - 4)*(x5 + 16/127)^2, -(2*x2 - 4)*(x5 + 16/127)^2, 1/10 - ((2*x5 + 32/127)*((x1 - 2)^2 + (x2 - 2)^2 - 1/16))/10 - (2*x5 + 32/127)*(x3*(2*x1 - 4) + x4*(2*x2 - 4))];
        fx_z=fx_z(x1,x2,x3,x4);

        %         fu_z=double(subs(fu_z,[x1,x2,x3,x4,x5,u1,u2],[x',u']));
        fu_z=[0,0];
        fx=eye(4)+dt*[0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
        fu=dt*[0 0;0 0;1 0;0 1];
        fbarx=[fx zeros(4,1);fx_z];
        fbaru=[fu;fu_z];
    else
        [fbarx, fbaru]=deal([]);
    end

end
