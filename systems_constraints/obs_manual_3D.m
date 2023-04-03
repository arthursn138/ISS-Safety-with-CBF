function [h,obs_info]=obs_manual_3D()
    [sphere_x,sphere_y,sphere_z]=sphere;
    x_ax(1)=2;  y_ax(1)=2.2;  z_ax(1)=1;   r(1)=1;
    x_ax(2)=-2; y_ax(2)=-2.2; z_ax(2)=1;   r(2)=1;
    x_ax(3)=4;  y_ax(3)=1;    z_ax(3)=-1;  r(3)=1;
    x_ax(4)=0;  y_ax(4)=3;    z_ax(4)=3.5; r(4)=1;
    x_ax(5)=-3; y_ax(5)=2;    z_ax(5)=-2;  r(5)=1;
    x_ax(6)=6;  y_ax(6)=0;    z_ax(6)=-4;  r(6)=1;
    x_ax(7)=1;  y_ax(7)=-3;   z_ax(7)=-3;  r(7)=1;
    x_ax(8)=1;  y_ax(8)=-3;   z_ax(8)=4;   r(8)=1;
    x_ax(9)=1;  y_ax(9)=-3;   z_ax(9)=0;  r(9)=1;
    
%     x_ax(10)=12;  y_ax(10)=2.2;  z_ax(10)=1;   r(10)=1;
%     x_ax(11)=8; y_ax(11)=-2.2; z_ax(11)=1;   r(11)=1;
%     x_ax(12)=14;  y_ax(12)=1;    z_ax(12)=-1;  r(12)=1;
%     x_ax(13)=10;  y_ax(13)=3;    z_ax(13)=3.5; r(13)=1;
%     x_ax(14)=7; y_ax(14)=2;    z_ax(14)=-2;  r(14)=1;
%     x_ax(15)=16;  y_ax(15)=0;    z_ax(15)=-4;  r(15)=1;
%     x_ax(16)=11;  y_ax(16)=-3;   z_ax(16)=-3;  r(16)=1;
%     x_ax(17)=11;  y_ax(17)=-3;   z_ax(17)=4;   r(17)=1;
%     x_ax(18)=11;  y_ax(18)=-3;   z_ax(18)=0;   r(18)=1;
%     x_ax(19)=5;  y_ax(19)=-0.5;   z_ax(19)=5;   r(19)=1;
    
    obs_info=[x_ax',y_ax',z_ax',r'];
    
    num_obs=size(obs_info,1);
    
    h=cell(num_obs,1);
    for ii=1:num_obs
        h{ii}= @(x) (x(10)-x_ax(ii))^2+ (x(11)-y_ax(ii))^2+(x(12)-z_ax(ii))^2-r(ii)^2;
    end
%     
%     figure(10000)
%     for ii=1:num_obs
%     spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
%     hold on;
%     spheres.EdgeColor = 'k';
%     spheres.FaceColor = '#A2142F';
%     spheres.LineStyle = ':';
%     spheres.FaceLighting = 'flat';
%     spheres.FaceAlpha= 1;
%     end
%     axis equal;   
%     xlabel('x'); ylabel('y'); zlabel('z');

    
end
    
