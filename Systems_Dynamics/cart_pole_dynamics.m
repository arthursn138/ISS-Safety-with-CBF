function [f,fx,fu,fxx,fxu,fuu,x,u] = cart_pole_dynamics(dt,deriv_bool)
    nx = 4; nu = 1; % state and input dimentions
    x = sym('x',[nx 1]);
    u = sym('u');
    nx=length(x);
    
    % systems paramters
    mp = 0.05;      % pole mass (kg)
    mc = 1;         % cart mass (kg)
    g = 9.8;        % gravitational force (m/s^2)
    l = 2;        % pole length (m)

%     xc = x(1); theta = x(2); xc_dot = x(3); theta_dot= x(4);
%     xc_dot = xc_dot;
%     theta_dot = theta_dot;
%     xc_ddot = @(mp,mc,g,l) (mp*sin(theta)*(l*theta_dot^2+g*cos(theta))+u)/(mc+mp*(sin(theta))^2);
%     theta_ddot = @(mp,mc,g,l)(-mp*l*theta_dot^2*cos(theta)*sin(theta)-(mc+mp)*g*sin(theta)-cos(theta)*u)/(l*(mc+mp*(sin(theta))^2));

%     f_cont_sym = @(x,u) [xc_dot;theta_dot;xc_ddot;theta_ddot];
    
    f = @(x,u) x(1:nx) + dt * [x(3);
                               x(4);
                               (mp*sin(x(2))*(l*x(4)^2+g*cos(x(2)))+u)/(mc+mp*(sin(x(2)))^2);
                               (-mp*l*x(4)^2*cos(x(2))*sin(x(2))-(mc+mp)*g*sin(x(2))-cos(x(2))*u)/(l*(mc+mp*(sin(x(2)))^2))];
    
    [fx, fu, fxx, fxu, fuu]=deal([]);

    if deriv_bool
        fx = @(x,u) eye(nx)+dt*[ 0, 0, 1, 0;
                                 0, 0, 0, 1;
                                 0, (mp*cos(x(2))*(l*x(4)^2 + g*cos(x(2))) - g*mp*sin(x(2))^2)/(mp*sin(x(2))^2 + mc) - (2*mp*cos(x(2))*sin(x(2))*(u + mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2)))))/(mp*sin(x(2))^2 + mc)^2, 0, (2*l*mp*x(4)*sin(x(2)))/(mp*sin(x(2))^2 + mc);
                                 0, (- l*mp*x(4)^2*cos(x(2))^2 + l*mp*x(4)^2*sin(x(2))^2 - g*(mc + mp)*cos(x(2)) + u*sin(x(2)))/(l*(mp*sin(x(2))^2 + mc)) + (2*mp*cos(x(2))*sin(x(2))*(l*mp*cos(x(2))*sin(x(2))*x(4)^2 + u*cos(x(2)) + g*sin(x(2))*(mc + mp)))/(l*(mp*sin(x(2))^2 + mc)^2), 0, -(2*mp*x(4)*cos(x(2))*sin(x(2)))/(mp*sin(x(2))^2 + mc)];
         fu = @(x,u) [0;0;1/(mp*sin(x(2))^2 + mc);-cos(x(2))/(l*(mp*sin(x(2))^2 + mc))]*dt;
        fx_temp = fx(x,u);
        for ii=1:nx
                fxx_temp = jacobian(fx_temp(ii,:), x);
                fxx{ii}=matlabFunction(fxx_temp,'Vars',{x,u});
                fxu_temp = jacobian(fx_temp(ii,:), u).';
                fxu{ii}=matlabFunction(fxu_temp,'Vars',{x,u});
        end
        fu_temp = fu(x,u);
        for ii=1:nu
                fuu_temp = jacobian(fu_temp(ii,:), u);
                fuu{ii}=matlabFunction(fuu_temp,'Vars',{x,u});
        end
%          fxx{1} = @(x,u) dt * [0, 0, 0, 0;0, 0, 0, 0;0, 0, 0, 0;0, 0, 0, 0];
%          fxx{2} = @(x,u) dt * [0, 0, 0, 0;0, 0, 0, 0;0, 0, 0, 0;0, 0, 0, 0];
%          fxx{3} = @(x,u) dt * [0,0,0,0;
%                                0,(2*mp*sin(x(2))^2*(u + mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2)))))/(mp*sin(x(2))^2 + mc)^2 - (2*mp*cos(x(2))^2*(u + mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2)))))/(mp*sin(x(2))^2 + mc)^2 - (mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2))) + 3*g*mp*cos(x(2))*sin(x(2)))/(mp*sin(x(2))^2 + mc) - (4*mp*cos(x(2))*sin(x(2))*(mp*cos(x(2))*(l*x(4)^2 + g*cos(x(2))) - g*mp*sin(x(2))^2))/(mp*sin(x(2))^2 + mc)^2 + (8*mp^2*cos(x(2))^2*sin(x(2))^2*(u + mp*sin(x(2))*(l*x(4)^2 + g*cos(x(2)))))/(mp*sin(x(2))^2 + mc)^3, 0, (2*l*mp*x(4)*cos(x(2)))/(mp*sin(x(2))^2 + mc) - (4*l*mp^2*x(4)*cos(x(2))*sin(x(2))^2)/(mp*sin(x(2))^2 + mc)^2;
%                                0,0, 0,0;
%                                0,(2*l*mp*x(4)*cos(x(2)))/(mp*sin(x(2))^2 + mc) - (4*l*mp^2*x(4)*cos(x(2))*sin(x(2))^2)/(mp*sin(x(2))^2 + mc)^2, 0,(2*l*mp*sin(x(2)))/(mp*sin(x(2))^2 + mc)];
%          fxx{4} = @(x,u) dt *[0,0,0,0;
%                               0,(4*l*mp*cos(x(2))*sin(x(2))*x(4)^2 + u*cos(x(2)) + g*sin(x(2))*(mc + mp))/(l*(mp*sin(x(2))^2 + mc)) + (2*mp*cos(x(2))^2*(l*mp*cos(x(2))*sin(x(2))*x(4)^2 + u*cos(x(2)) + g*sin(x(2))*(mc + mp)))/(l*(mp*sin(x(2))^2 + mc)^2) - (2*mp*sin(x(2))^2*(l*mp*cos(x(2))*sin(x(2))*x(4)^2 + u*cos(x(2)) + g*sin(x(2))*(mc + mp)))/(l*(mp*sin(x(2))^2 + mc)^2) - (4*mp*cos(x(2))*sin(x(2))*(- l*mp*x(4)^2*cos(x(2))^2 + l*mp*x(4)^2*sin(x(2))^2 - g*(mc + mp)*cos(x(2)) + u*sin(x(2))))/(l*(mp*sin(x(2))^2 + mc)^2) - (8*mp^2*cos(x(2))^2*sin(x(2))^2*(l*mp*cos(x(2))*sin(x(2))*x(4)^2 + u*cos(x(2)) + g*sin(x(2))*(mc + mp)))/(l*(mp*sin(x(2))^2 + mc)^3), 0, (2*mp*x(4)*sin(x(2))^2)/(mp*sin(x(2))^2 + mc) - (2*mp*x(4)*cos(x(2))^2)/(mp*sin(x(2))^2 + mc) + (4*mp^2*x(4)*cos(x(2))^2*sin(x(2))^2)/(mp*sin(x(2))^2 + mc)^2;
%                               0,0,0,0;
%                               0,(2*mp*x(4)*sin(x(2))^2)/(mp*sin(x(2))^2 + mc) - (2*mp*x(4)*cos(x(2))^2)/(mp*sin(x(2))^2 + mc) + (4*mp^2*x(4)*cos(x(2))^2*sin(x(2))^2)/(mp*sin(x(2))^2 + mc)^2, 0, -(2*mp*cos(x(2))*sin(x(2)))/(mp*sin(x(2))^2 + mc)];
% 
%          %fux transposed:
%          fxu= @(x,u) dt * [0, 0, 0, 0;
%                            0, 0, 0, 0;
%                            0,-(2*mp*cos(x(2))*sin(x(2)))/(mp*sin(x(2))^2 + mc)^2, 0, 0;
%                            0, sin(x(2))/(l*(mp*sin(x(2))^2 + mc)) + (2*mp*cos(x(2))^2*sin(x(2)))/(l*(mp*sin(x(2))^2 + mc)^2), 0, 0];
% 
%          fuu = @(x,u) zeros(4,1);
    end
end