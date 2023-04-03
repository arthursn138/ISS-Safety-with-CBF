function [f,fx,fu,fxx,fxu,fuu,x,u] = Diff_wheel_dynamics(dt,deriv_bool)
%     syms x1 x2 x3 u1 u2
%     syms x u
%     x=[x(1);x(2);x(3)];
%     u=[u(1);u(2)];
    n = 3; m = 2; % state and input dimentions
    x = sym('x',[n 1]);
    u = sym('u',[m 1]);
    n=length(x);
    d = 0.2;        %tread of the robot
    r = 0.2;       %wheel rad of the right wheel
%     rleft = 0.2;       %wheel rad of the left wheel
    v = (r*u(1)+r*u(2))/2;
    omega = (r*u(1)-r*u(2))/(2*d);
    
    xdot = v*cos(x(3));
    ydot = v*sin(x(3));
    thetadot = omega;

%     f_cont = [xdot;ydot;thetadot];

%     f = @(x,u) [x(1);x(2);x(3)] + dt * [xdot;ydot;thetadot];
    
    f = @(x,u) x(1:n) + dt * [(r*u(1)+r*u(2))/2*cos(x(3));(r*u(1)+r*u(2))/2*sin(x(3));(r*u(1)-r*u(2))/(2*d)];
    
    [fx, fu, fxx, fxu, fuu]=deal([]);
    if deriv_bool
        fx = @(x,u) eye(n)+dt*[0 0 -(r*u(1)+r*u(2))/2*sin(x(3)); 0 0 (r*u(1)+r*u(2))/2*cos(x(3)); 0 0 0];
        fu = @(x,u) [r*cos(x(3))/2,r*cos(x(3))/2;r*sin(x(3))/2,r*sin(x(3))/2;...
        r/(2*d) -r/(2*d)]*dt;
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

