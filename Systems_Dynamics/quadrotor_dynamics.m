function [f,fx,fu,x,u] = quadrotor_dynamics(dt,deriv_bool)
%     syms x(1).. x(12) u(1) .. u(4)
%     syms x u
%     x=[x(1);...;x(12)];
%     u=[u(1);...;u(4)];
    n = 12; m = 4;% state and input dimentions
    x = sym('x',[n 1]);
    u = sym('u',[m 1]);
    n=length(x);

    g = 9.81;  % gravity constant
    mass =1;   % quad mass
    Ix =1;  % quad x-axis inertia
    Iy =1;  % quad y-axis inertia
    Iz =1;  % quad z-axis inertia
    tau_wx = 0;% torque due to wind on the x-axis
    tau_wy = 0;% torque due to wind on the y-axis
    tau_wz = 0;% torque due to wind on the z-axis
    f_wx = 0;  % force due to wind on the x-axis
    f_wy = 0;  % force due to wind on the y-axis
    f_wz = 0;  % force due to wind on the z-axis
    
%     syms g mass Ix Iy Iz tau_wx tau_wy tau_wz f_wx w_wy f_wz
    phi = x(1);% roll angle (x-axis) in earth frame
    theta = x(2);% pitch angle (y-axis) in earth frame
    psi = x(3);% yaw angle (z-axis) in earth frame
    phi_rate = x(4);   % roll rate/angular velocity in the body fram
    theta_rate = x(5); % pitch rate/angular velocity in the body fram
    psi_rate = x(6);   % yaw rate/angular velocity in the body fram
    v_x = x(7);% x-linear velocity in the body frame
    v_y = x(8);% y-linear velocity in the body frame
    v_z = x(9);% z-linear velocity in the body frame
    x_axis = x(10);    % x position in the earth frame
    y_axis = x(11);    % y position in the earth frame
    z_axis = x(12);    % z position in the earth frame
    ft = u(1);% thrust in the body frame
    tau_x = u(2);% torque along x in the body frame
    tau_y = u(3);% torque along y in the body frame
    tau_z = u(4);% torque along z in the body frame
    
    phi_dot = phi_rate+psi_rate*cos(phi)*tan(theta)+theta_rate*sin(phi)*tan(theta);
    theta_dot = theta_rate*cos(phi)-psi_rate*sin(phi);
    psi_dot = psi_rate*cos(phi)/cos(theta)+theta_rate*sin(phi)/cos(theta);
    phi_rate_dot = (Iy-Iz)/Ix*theta_rate*psi_rate+(tau_x+tau_wx)/Ix;
    theta_rate_dot = (Iz-Ix)/Iy*phi_rate*psi_rate+(tau_y+tau_wy)/Iy;
    psi_rate_dot = (Ix-Iy)/Iz*phi_rate*theta_rate+(tau_z+tau_wz)/Iz;
    v_x_dot = psi_rate*v_y-theta_rate*v_z-g*sin(theta)+f_wx/mass;
    v_y_dot = phi_rate*v_z-psi_rate*v_x+g*sin(phi)*cos(theta)+f_wy/mass;
    v_z_dot = theta_rate*v_x-phi_rate*v_y+g*cos(theta)*cos(phi)+(f_wz-ft)/mass;
    x_axis_dot = v_z*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) - v_y*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta)) + v_x*(cos(psi)*cos(theta));
    y_axis_dot = v_y*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) - v_z*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta)) + v_x*(cos(theta)*sin(psi));
    z_axis_dot = v_z*cos(phi)*cos(theta)-v_x*sin(theta)+v_y*cos(theta)*sin(phi);
% symbolic cont. time dynamics for generating fx nd fu first time only 
%     f_cont_sym = @(x,u) [phi_dot;theta_dot;psi_dot;phi_rate_dot;theta_rate_dot;psi_rate_dot;v_x_dot;v_y_dot;v_z_dot;x_axis_dot;y_axis_dot;z_axis_dot];
    
    f = @(x,u) x(1:n) + dt * [x(4) + x(6)*cos(x(1))*tan(x(2)) + x(5)*sin(x(1))*tan(x(2));
                              x(5)*cos(x(1)) - x(6)*sin(x(1));
                              (x(6)*cos(x(1)))/cos(x(2)) + (x(5)*sin(x(1)))/cos(x(2));
                              (tau_wx + u(2))/Ix + (x(5)*x(6)*(Iy - Iz))/Ix;
                              (tau_wy + u(3))/Iy - (x(4)*x(6)*(Ix - Iz))/Iy;
                              (tau_wz + u(4))/Iz + (x(4)*x(5)*(Ix - Iy))/Iz;
                              x(6)*x(8) - x(5)*x(9) + f_wx/mass - g*sin(x(2));
                              x(4)*x(9) - x(6)*x(7) + f_wy/mass + g*cos(x(2))*sin(x(1));
                              x(5)*x(7) - x(4)*x(8) - (u(1) - f_wz)/mass + g*cos(x(1))*cos(x(2));
                              x(9)*(sin(x(1))*sin(x(3)) + cos(x(1))*cos(x(3))*sin(x(2))) - x(8)*(cos(x(1))*sin(x(3)) - cos(x(3))*sin(x(1))*sin(x(2))) + x(7)*cos(x(2))*cos(x(3));
                              x(8)*(cos(x(1))*cos(x(3)) + sin(x(1))*sin(x(2))*sin(x(3))) - x(9)*(cos(x(3))*sin(x(1)) - cos(x(1))*sin(x(2))*sin(x(3))) + x(7)*cos(x(2))*sin(x(3));
                              x(9)*cos(x(1))*cos(x(2)) - x(7)*sin(x(2)) + x(8)*cos(x(2))*sin(x(1))];
    
    if deriv_bool
% fx = eye(n)+dt*[0 0 -v*sin(x(3)); 0 0 v*cos(x(3)); 0 0 0];
  fx = @(x,u) eye(n)+dt* [x(5)*cos(x(1))*tan(x(2)) - x(6)*sin(x(1))*tan(x(2)), x(6)*cos(x(1))*(tan(x(2))^2 + 1) + x(5)*sin(x(1))*(tan(x(2))^2 + 1), 0, 1, sin(x(1))*tan(x(2)), cos(x(1))*tan(x(2)), 0, 0, 0, 0, 0, 0;
                        -x(6)*cos(x(1))-x(5)*sin(x(1)), 0, 0, 0, cos(x(1)), -sin(x(1)), 0, 0, 0, 0, 0, 0;
                        (x(5)*cos(x(1)))/cos(x(2))-(x(6)*sin(x(1)))/cos(x(2)), (x(6)*cos(x(1))*sin(x(2)))/cos(x(2))^2+(x(5)*sin(x(1))*sin(x(2)))/cos(x(2))^2, 0, 0, sin(x(1))/cos(x(2)), cos(x(1))/cos(x(2)), 0, 0, 0, 0, 0, 0;
                        0, 0, 0, 0, (x(6)*(Iy-Iz))/Ix, (x(5)*(Iy-Iz))/Ix, 0, 0, 0, 0, 0, 0;
                        0, 0, 0, -(x(6)*(Ix-Iz))/Iy, 0, -(x(4)*(Ix-Iz))/Iy, 0, 0, 0, 0, 0, 0;
                        0, 0, 0, (x(5)*(Ix-Iy))/Iz, (x(4)*(Ix-Iy))/Iz, 0, 0, 0, 0, 0, 0, 0;
                        0, -g*cos(x(2)), 0, 0, -x(9), x(8), 0, x(6), -x(5), 0, 0, 0;
                        g*cos(x(1))*cos(x(2)), -g*sin(x(1))*sin(x(2)), 0, x(9), 0, -x(7), -x(6), 0, x(4), 0, 0, 0;
                        -g*cos(x(2))*sin(x(1)), -g*cos(x(1))*sin(x(2)), 0, -x(8), x(7), 0, x(5), -x(4), 0, 0, 0, 0;
                        x(8)*(sin(x(1))*sin(x(3))+cos(x(1))*cos(x(3))*sin(x(2)))+x(9)*(cos(x(1))*sin(x(3))-cos(x(3))*sin(x(1))*sin(x(2))), x(9)*cos(x(1))*cos(x(2))*cos(x(3))-x(7)*cos(x(3))*sin(x(2))+x(8)*cos(x(2))*cos(x(3))*sin(x(1)), x(9)*(cos(x(3))*sin(x(1))-cos(x(1))*sin(x(2))*sin(x(3)))-x(8)*(cos(x(1))*cos(x(3)) + sin(x(1))*sin(x(2))*sin(x(3))) - x(7)*cos(x(2))*sin(x(3)), 0,  0, 0, cos(x(2))*cos(x(3)), cos(x(3))*sin(x(1))*sin(x(2)) - cos(x(1))*sin(x(3)), sin(x(1))*sin(x(3)) + cos(x(1))*cos(x(3))*sin(x(2)), 0, 0, 0;
                        -x(8)*(cos(x(3))*sin(x(1)) - cos(x(1))*sin(x(2))*sin(x(3))) - x(9)*(cos(x(1))*cos(x(3)) + sin(x(1))*sin(x(2))*sin(x(3))), x(9)*cos(x(1))*cos(x(2))*sin(x(3)) - x(7)*sin(x(2))*sin(x(3)) + x(8)*cos(x(2))*sin(x(1))*sin(x(3)), x(9)*(sin(x(1))*sin(x(3)) + cos(x(1))*cos(x(3))*sin(x(2))) - x(8)*(cos(x(1))*sin(x(3)) - cos(x(3))*sin(x(1))*sin(x(2))) + x(7)*cos(x(2))*cos(x(3)), 0,  0, 0, cos(x(2))*sin(x(3)), cos(x(1))*cos(x(3)) + sin(x(1))*sin(x(2))*sin(x(3)), cos(x(1))*sin(x(2))*sin(x(3)) - cos(x(3))*sin(x(1)), 0, 0, 0;
                        x(8)*cos(x(1))*cos(x(2)) - x(9)*cos(x(2))*sin(x(1)),   - x(7)*cos(x(2)) - x(9)*cos(x(1))*sin(x(2)) - x(8)*sin(x(1))*sin(x(2)), 0, 0,  0, 0, -sin(x(2)), cos(x(2))*sin(x(1)), cos(x(1))*cos(x(2)), 0, 0, 0];
 
 fu = @(x,u) dt*[ 0,    0,    0,    0;
                  0,    0,    0,    0;
                  0,    0,    0,    0;
                  0, 1/Ix,    0,    0;
                  0,    0, 1/Iy,    0;
                  0,    0,    0, 1/Iz;
                  0,    0,    0,    0;
                  0,    0,    0,    0;
            -1/mass,    0,    0,    0;
                  0,    0,    0,    0;
                  0,    0,    0,    0;
                  0,    0,    0,    0];
    else
 [fx, fu]=deal([]);
    end
end

