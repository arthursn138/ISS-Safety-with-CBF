function [f,fx,fu,x,u] = toy_nonaffine_dynamics(dt,deriv_bool)
%     syms x1 x2 x3 u1 u2
%     syms x u
%     x=[x(1);x(2);x(3)];
%     u=[u(1);u(2)];
    n = 3; m = 2; % state and input dimentions
    x = sym('x',[n 1]);
    u = sym('u',[m 1]);
    n=length(x);

%     xdot = u(1)*cos(x(3))+u(1)*u(2)/2;
%     ydot = u(2)^2*sin(x(3));
%     zdot = u(1)+u(2);

%     f_cont = @(x,u)[xdot;ydot;zdot];
    
    f = @(x,u) x(1:n) + dt * [u(1)*cos(x(3))+u(1)*u(2)/2;u(2)^2*sin(x(3));u(1)+u(2)];
    
    if deriv_bool
         fx = @(x,u) eye(n)+dt*[0, 0, -u(1)*sin(x(3)); 0, 0, u(2)^2*cos(x(3)); 0 0 0];
        fu = @(x,u) [u(2)/2 + cos(x(3)), u(1)/2;0, 2*u(2)*sin(x(3));1, 1]*dt;
    else
        [fx, fu]=deal([]);
    end
end

