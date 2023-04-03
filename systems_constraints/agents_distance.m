function [h]=agents_distance()

delta =0.5;

h{1} = @(x) (x(1)-x(5))^2 + (x(2) - x(6))^2 - delta^2;

%     for ii=1:num_agents
%         h{ii}= @(x) (x(1)-x_ax(ii))^2+ (x(2)-y_ax(ii))^2-r(ii)^2;
%     end
    

end
    
