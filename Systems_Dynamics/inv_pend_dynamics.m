function [f,fx,fu,x,u] = inv_pend_dynamics(dt,deriv_bool)
    nx = 2; nu = 1; % state and input dimentions
    x = sym('x',[nx 1]);
    u = sym('u');
    uncertinity_perc=00/100; c=1-uncertinity_perc;
    % systems paramters
    m = 1.5*c;      % mass (kg)
    b = 0.15*c;    % friction (kg)
    g = 9.81;   % gravitational force (m/s^2)
    l = 0.75*c;    % length (m)
    I=m*l^2;    % inertia


%Eq of motion: I*theta_ddot + b * theta_dot - m*g*l * sign(theta) = u; 

%     theta = x(1); theta_dot = x(2);
%     theta_dot = theta_dot;
%     theta_ddot= @(m,b,g,l,I) (-b * theta_dot + m*g*l * sin(theta) + u)/I;
%     f_cont_sym = @(x,u) [theta_dot;theta_ddot];
    
    f = @(x,u) x(1:nx) + dt * [x(2);
                               (-b * x(2) + m*g*l * sin(x(1)) + u)/I];
    
    if deriv_bool
        fx = @(x,u) eye(nx)+dt*[ 0, 1;
                                 (g*l*m*cos(x(1)))/I, -b/I];
        fu = @(x,u) [0;1/I]*dt;
    else
        [fx, fu]=deal([]);
    end
end
