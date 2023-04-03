function [h,w_max, gamma_gs,theta_max, delta_max, Tmin, Tmax]=rocket_const(xf)

mdry=0.0;
Tmin=0.5; Tmax=5;% Tmax = 5;
delta_max=10;
theta_max=20;
gamma_gs=10;
H23=[0 0;1 0; 0 1];
e1=[1;0;0];
w_max=30;%30*pi/180;
% x(1)>= mdry
% r_I*e1 - tan(gamma_gs*norm(H23*r_I)) >= 0
% cos(theta_max) <= 1 - 2*(q2^2 + q3^2)
% norm(w_B) <= w_max
rel=0;
%  h{1} = @(x) x(1) + 1e-1; h{2} = @(x) x(2) + 0.1; h{3} = @(x) x(3)+ 1e-1;
% h{1} = @(x) rel+ x(1) - mdry;

h{1} = @(x) w_max^2 - norm(x(11:13))^2; % fine

% h{1} = @(x) x(1) - tand(gamma_gs)*norm([x(2);x(3)]); % fine
h{2} = @(x) 1 - 2 * (x(9)^2 + x(10)^2) - cosd(theta_max);

% h{3} = @(x) Tmax - norm(x(14:16));
% h{4} = @(x) (x(14)^2+x(15)^2+x(16)^2) - Tmin^2;
% h{3} = @(x) 1e-2 + x(14)-(cosd(delta_max))^2*norm(x(14:16));


% h{2} = @(x) e1.'*[x(14);x(15);x(16)]-cosd(delta_max)*norm(x(14:16));
% h{2} = @(x) e1.'*[x(14);x(15);x(16)]-cosd(delta_max*norm(x(14:16)));



% h{3} = @(x) 1e-10 + (e1.'*x(1:3)) - tan(gamma_gs*norm(H23.'*x(1:3)));
% h{3} = @(x) rel+1e-10 + x(1:3).' * e1 - tan(gamma_gs)*norm(H23.'*x(1:3));
end
    
