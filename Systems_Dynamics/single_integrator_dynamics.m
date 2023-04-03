%% Dynamics of Planar Double Integrator
function [f,fx,fu,fxx,fxu,fuu,x,u]=single_integrator_dynamics(dt,deriv_bool)
    syms x1 x2 u1 u2
    x=[x1;x2];
    u=[u1;u2];
    n=2; m=2;
% x1dot= u1;
% x2dot= u2;
% xdot=f(x,u)
% f=[u(1); u(2)];
    f=@(x,u) [x(1);x(2)]+dt*[u(1); u(2)];

    [fx, fu, fxx, fxu, fuu]=deal([]);

    if deriv_bool
        % dynamics gradients
%         fx=[0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
        fx=@(x,u)  eye(n)+dt*[0 0;0 0];
%         fu=[0 0;0 0;1 0;0 1];
        fu=@(x,u)  dt*[1 0;0 1];
        fx_temp = fx(x,u);
        for ii=1:n
                fxx_temp = jacobian(fx_temp(ii,:), x);
                fxx{ii}=matlabFunction(fxx_temp,'Vars',{x,u});
                fxu_temp = jacobian(fx_temp(ii,:), u).';
                fxu{ii}=matlabFunction(fxu_temp,'Vars',{x,u});
        end
        fu_temp = fu(x,u);
        for ii=1:m
                fuu_temp = jacobian(fu_temp(ii,:), u);
                fuu{ii}=matlabFunction(fuu_temp,'Vars',{x,u});
        end
    end

end


