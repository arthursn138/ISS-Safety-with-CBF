% % % cart_position = X(1,:);
% % % theta = X(2,:);
% % % x_pole=[]; y_pole=[]; L =2;
% % % for ii=1:N
% % % x_pole(ii)=cart_position(ii)+L*cos(theta(ii)-pi/2);
% % % y_pole(ii)=L*sin(theta(ii)-pi/2);
% % % end
% % % 
% % % fig =figure(100); 
% % % axis tight manual % this ensures that getframe() returns a consistent size
% % % for idx = 1:1:N
% % %     plot(cart_position(idx),0,'s',[cart_position(idx),x_pole(idx)],[0;y_pole(idx)],'-',x_pole(1:idx),y_pole(1:idx),'--',x_pole(idx),y_pole(idx),'o','LineWidth',1.5); 
% % %     rectangle('Position',[cart_limit,-0.1,100,.2],'FaceColor',[0 .5 .5]);
% % %     rectangle('Position',[-100,-0.1,100-cart_limit,.2],'FaceColor',[0 .5 .5]);
% % %     axis([-3 3 -3 3]);
% % %     drawnow
% % %     frame = getframe(fig);
% % %     im{idx} = frame2im(frame);
% % % end

fig =figure(100); 
axis tight manual % this ensures that getframe() returns a consistent size
kk=1;
for ll=1:num_trials

% % %     if (ll==3)
% % %         continue;
% % %     end
fig =figure(100); 
axis square;
X1_bas(ll,:);
    obstacles=ALLOBS(1+num_obs*(ll-1):ll*num_obs,:);
subplot(1,2,1);
for ii=1:size(obstacles,1)
   rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0 0 0]); hold on
end
scatter(Xf(1,ll),Xf(2,ll),80,'filled','hg'); 
scatter(X0(1,ll),X0(2,ll),50,'filled','or');
plot(X1_bas(ll,:),X2_bas(ll,:),'LineWidth',0.5); 
axis ([-3.5 3.5 -3.5 3.5]); axis square;
title('DBaS-DDP','FontName','Times New Roman','Interpreter','latex');
subplot(1,2,2);
for ii=1:size(obstacles,1)
   rectangle('Position',obstacles(ii,:),'Curvature',[1 1],'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0 0 0]); hold on
end
scatter(Xf(1,ll),Xf(2,ll),80,'filled','hg'); 
scatter(X0(1,ll),X0(2,ll),50,'filled','or');
plot(X1_pen(ll,:),X2_pen(ll,:),'LineWidth',0.5); 
xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
axis ([-3.5 3.5 -3.5 3.5]); axis square;
title('Penalty-DDP','FontName','Times New Roman','Interpreter','latex');
drawnow
frame = getframe(fig);
im{kk} = frame2im(frame);
close all
    kk=kk+1;
end
filename = '/home/halmubarak3/Desktop/diff_wheel_animate10.gif'; % Specify the output file name
for idx = 1:num_trials
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.75);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.75);
    end
end
