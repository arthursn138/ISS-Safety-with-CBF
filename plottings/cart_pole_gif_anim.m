cart_position = X(1,:);
theta = X(2,:);
x_pole=[]; y_pole=[]; L =2;
for ii=1:N
x_pole(ii)=cart_position(ii)+L*cos(theta(ii)-pi/2);
y_pole(ii)=L*sin(theta(ii)-pi/2);
end
U(150)=U(149)
fig =figure(100); 
axis tight manual % this ensures that getframe() returns a consistent size
fig.Position = [200  100  1300  500];
t=dt;
for idx = 1:1:N
    subplot(1,2,1)
    plot(cart_position(idx),0,'s',[cart_position(idx),x_pole(idx)],[0;y_pole(idx)],'-',x_pole(1:idx),y_pole(1:idx),'--',x_pole(idx),y_pole(idx),'o','LineWidth',1.5); 
    rectangle('Position',[cart_limit+0.05,-0.1,100,.2],'FaceColor',[0 .5 .5]);
    rectangle('Position',[-100,-0.1,100-cart_limit-0.05,.2],'FaceColor',[0 .5 .5]);
    axis([-3 3 -3 3]);
    title('DBaS DDP','FontName','Times New Roman','Interpreter','latex');
        text(1.75,2,'t='); text(2,2,string(t)); text(2.4,2,'sec')
    subplot(1,2,2)
    plot(T(1:idx),U(1:idx),'-','LineWidth',1.5); 
    axis([0 3 -70 70]);
    title('DBaS DDP Control','FontName','Times New Roman','Interpreter','latex');        
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
    t=t+dt;
end
%%
filename = 'C:\Users\HasMBK\Desktop/cartpole_dbasddp_animate.gif'; % Specify the output file name
ppause=0.05;
for idx = 1:N
    [A,map] = rgb2ind(im{idx},256);
    if idx==N
        ppause=1;
    end
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',ppause);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',ppause);
    end
end