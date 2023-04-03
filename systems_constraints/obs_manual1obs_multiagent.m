function [h,obs_loc,circ]=obs_manual1obs_multiagent()
    num_obs=8;

    x_ax(1)=0; y_ax(1)=2; r(1)=1; 
    circ(1,:) = [x_ax(1)-r(1) y_ax(1)-r(1) 2*r(1) 2*r(1)];
    
    x_ax(2)=0; y_ax(2)=4; r(2)=1; 
    circ(2,:) = [x_ax(2)-r(2) y_ax(2)-r(2) 2*r(2) 2*r(2)];
    
    x_ax(3)=0; y_ax(3)=6; r(3)=1; 
    circ(3,:) = [x_ax(3)-r(3) y_ax(3)-r(3) 2*r(3) 2*r(3)];
    
    x_ax(4)=0; y_ax(4)=8; r(4)=1; 
    circ(4,:) = [x_ax(4)-r(4) y_ax(4)-r(4) 2*r(4) 2*r(4)];

    x_ax(5)=0; y_ax(5)=-2; r(5)=1; 
    circ(5,:) = [x_ax(5)-r(5) y_ax(5)-r(5) 2*r(5) 2*r(5)];

    x_ax(6)=0; y_ax(6)=-4; r(6)=1; 
    circ(6,:) = [x_ax(6)-r(6) y_ax(6)-r(6) 2*r(6) 2*r(6)];

    x_ax(7)=0; y_ax(7)=-6; r(7)=1; 
    circ(7,:) = [x_ax(7)-r(7) y_ax(7)-r(7) 2*r(7) 2*r(7)];

    x_ax(8)=0; y_ax(8)=-8; r(8)=1; 
    circ(8,:) = [x_ax(8)-r(8) y_ax(8)-r(8) 2*r(8) 2*r(8)];

    
    obs_loc=[x_ax,y_ax,r];

 
    h=cell(num_obs,1);
    ind=1;
    for ii=1:num_obs
        h{ind}= @(x) (x(1)-x_ax(ii))^2+ (x(2)-y_ax(ii))^2-r(ii)^2;
        h{ind+1}= @(x) (x(5)-x_ax(ii))^2+ (x(6)-y_ax(ii))^2-r(ii)^2;
        ind=ind+2;
    end
    
%     figure(10000)
%     for ii=1:num_obs
%     rectangle('Position',circ(ii,:),'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
%     end
   
end
    
