% cart_position = X(1,:);
% theta = X(2,:);
% x_pole=[]; y_pole=[]; L =2;
% for ii=1:N
% x_pole(ii)=cart_position(ii)+L*cos(theta(ii)-pi/2);
% y_pole(ii)=L*sin(theta(ii)-pi/2);
% end
% 
% fig =figure(100); 
% axis tight manual % this ensures that getframe() returns a consistent size
% for idx = 1:1:N
%     plot(cart_position(idx),0,'s',[cart_position(idx),x_pole(idx)],[0;y_pole(idx)],'-',x_pole(1:idx),y_pole(1:idx),'--',x_pole(idx),y_pole(idx),'o','LineWidth',1.5); 
%     rectangle('Position',[cart_limit,-0.1,100,.2],'FaceColor',[0 .5 .5]);
%     rectangle('Position',[-100,-0.1,100-cart_limit,.2],'FaceColor',[0 .5 .5]);
%     axis([-3 3 -3 3]);
%     drawnow
%     frame = getframe(fig);
%     im{idx} = frame2im(frame);
% end

fig =figure(100); 
% axis tight manual % this ensures that getframe() returns a consistent size
% scatter(xf(1),xf(2),80,'filled','hg'); hold on;
% scatter(x0(1),x0(2),50,'filled','or');

% figure(100)
% for ii=1:size(obstacles,1)
%    rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); 
% end

% hold off;
k=0;
for ll=1:5:N
    k=k+1;
plot(X(1,1:ll),X(2,1:ll),'-',x0(1),x0(2),'o',X(1,ll),X(2,ll),'o',xf(1),xf(2),'gh','LineWidth',1.5); 
% plot(x0(1),x0(2),'ro',xf(1),xf(2),'gh','LineWidth',1.5); 
for ii=1:size(obstacles,1)
   rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840]); 
end
% scatter(x0(1),x0(2),50,'filled','or');
% plot(X(1,ll),X(2,ll),'o','Color','#0072BD','LineWidth',1);
title('Differential Wheeled Robot with DBaS-DDP','FontName','Times New Roman','Interpreter','latex');
axis ([-.5 5 -.7 5]); axis square;
drawnow
frame = getframe(fig);
im{k} = frame2im(frame);
end
filename = '/home/halmubarak3/Desktop/diff_wheel_manyspheres_dbasddp.gif'; % Specify the output file name
ppause=0.05;
for idx = 1:k
    [A,map] = rgb2ind(im{idx},256);
    if idx==k
        ppause=1;
    end
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',ppause);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',ppause);
    end
end
