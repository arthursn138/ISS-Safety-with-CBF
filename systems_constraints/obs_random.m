function [h,obs_loc,circ]=obs_random(num_obs)
%         close all
%     num_obs=10;
    
    for ii=1:num_obs
        x_ax(ii)=randn; 
        y_ax(ii)=randn; 
        r(ii)=rand; % abs(randn)/2
        circ(ii,:) = [x_ax(ii)-r(ii) y_ax(ii)-r(ii) 2*r(ii) 2*r(ii)];
    end
    obs_loc=[x_ax,y_ax,r];
        
    h=cell(num_obs,1);
    for ii=1:num_obs
        h{ii}= @(x) (x(1)-x_ax(ii))^2+ (x(2)-y_ax(ii))^2-r(ii)^2;
    end
    
%     figure(1)
%     for ii=1:num_obs
%     rectangle('Position',circ(ii,:),'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
%     end
   
end
    
