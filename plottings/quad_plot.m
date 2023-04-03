%% for quad - many sims 
function quad_plot(f_dyn,k_u,K_u,N,U,T,X,x0,xf,obs_info)
XX1=[]; XX2=[]; XX3=[]; BB=[];
% T=0:dt:dt*N-1*dt;
for kk=1:500
    [X_sim,U_sim]=simulate_ddp(f_dyn,k_u,K_u,x0,N,X,U);
% plotting (3D obstacle avoidance)
% figure(10000)
% plot3(X(10,:),X(11,:),X(12,:),'b-',x0(10),x0(11),x0(12),'bo',xf(10),xf(11),xf(12),'gx','LineWidth',1); hold on;
% plot3(X(10,:),X(11,:),X(12,:),'LineWidth',1); hold on;
% scatter3(xf(10),xf(11),xf(12),50,'filled','og','LineWidth',1); hold on
% scatter3(x0(10),x0(11),x0(12),50,'filled','or','LineWidth',1);
% plot3(X_sim(10,:),X_sim(11,:),X_sim(12,:),'-','Color','#D95319','LineWidth',.25); hold on;
% plot3(X_sim(10,:),X_sim(11,:),X_sim(12,:),'-','LineWidth',.25); hold on;
% grid on;
% 
% % 
%     num_obs=size(obs_info,1);
%     [sphere_x,sphere_y,sphere_z]=sphere;
%     x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
%     r = obs_info(:,4);
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
% xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
% ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
% zlabel('$z$','FontName','Times New Roman','Interpreter','latex');
% %     axis([-6 17 -6 6 -6 6])
% box on; axis square;
% box on; axis equal;

% figure(50)
% subplot(3,1,1);
% plot(0:dt:dt*N-2*dt,U_sim,'LineWidth',1.5); hold on
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$u$','FontName','Times New Roman','Interpreter','latex');
% subplot(3,1,2);
% plot(0:dt:dt*N-2*dt,V_sim,'LineWidth',1.5); hold on
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('$v$','FontName','Times New Roman','Interpreter','latex');
% subplot(3,1,3);
% plot(T,X_sim(13,:),'LineWidth',1.5); hold on;
% xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
% ylabel('Barrier State $w$','FontName','Times New Roman','Interpreter','latex');
XX1(kk,:)=X_sim(10,:);
XX2(kk,:)=X_sim(11,:);
XX3(kk,:)=X_sim(12,:);
BB(kk,:)=X_sim(13,:);
end

% Confidence region approx
figure(1000);
% for ii=1:size(XX1,1)
%     plot3(XX1(ii,:),XX2(ii,:),XX3(ii,:),'-','Color','#D95319','LineWidth',.25); hold on;
% end
    num_obs=size(obs_info,1);
    x_ax = obs_info(:,1); y_ax = obs_info(:,2); z_ax = obs_info(:,3);
    r = obs_info(:,4); [sphere_x,sphere_y,sphere_z]=sphere;
    for ii=1:num_obs
    spheres= surf(sphere_x*r(ii)+x_ax(ii), sphere_y*r(ii)+y_ax(ii), sphere_z*r(ii)+z_ax(ii));
    hold on;
    spheres.EdgeColor = 'k';
    spheres.FaceColor = '#A2142F';
    spheres.LineStyle = ':';
    spheres.FaceLighting = 'flat';
    spheres.FaceAlpha= 1;
    end
    axis equal;   
    
facecolor='#D95319'; alpha=0.01; edgecolor = 'none';
linestyle = 'none';
conf_perc=95/100; spc =1;
zstar=chi2inv(conf_perc,3);  % 3 is for 3D 
fig_num=1000;
plot_confidence_region(fig_num,XX1,XX2,XX3,N,zstar,spc,facecolor,alpha,edgecolor,linestyle)

plot3(X(10,:),X(11,:),X(12,:),'-','Color','#D95319','LineWidth',1.5); hold on;
plot3(x0(10),x0(11),x0(12),'ro','LineWidth',1)
plot3(xf(10),xf(11),xf(12),'gx','LineWidth',2); hold on;

% title('Safety Embedded Robust DDP','FontName','Times New Roman','Interpreter','latex');
xlabel('$x$','FontName','Times New Roman','Interpreter','latex');
ylabel('$y$','FontName','Times New Roman','Interpreter','latex');
zlabel('$z$','FontName','Times New Roman','Interpreter','latex');
% legend('Start','Target','Undisturbed Trajectory','$95\%$ Confidence Region','Unsafe Region','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)
grid on;

for ii=1:size(XX1,1)
plot3(XX1(ii,:),XX2(ii,:),XX3(ii,:),'-','Color',facecolor,'LineWidth',.25); hold on;
end
%
%% Confidence interval for barrier state:
BaS_mean = mean(BB);
BaS_std = std(BB);
zstar=sqrt(chi2inv(conf_perc,1));
% upper bound:
BaS_upper=BaS_mean + zstar*BaS_std;
BaS_lower=BaS_mean - zstar*BaS_std;
figure(51)
plot(T,X(13,:),'-','Color','#D95319','LineWidth',1.5); hold on;
% plot(T,BaS_mean,T,BaS_upper,T,BaS_lower); hold on;
conf_int = patch([T fliplr(T)],[min([BaS_upper;BaS_lower]) fliplr(max([BaS_upper;BaS_lower]))],'r')
conf_int.EdgeColor = edgecolor;
conf_int.LineStyle = linestyle;
conf_int.FaceColor = facecolor;
conf_int.FaceAlpha = 0.5;
xlabel('Time [sec]','FontName','Times New Roman','Interpreter','latex');
ylabel('Barrier State','FontName','Times New Roman','Interpreter','latex');
legend('Undisturbed Trajectory','$95\%$ Confidence Region','Unsafe Region','FontName','Times New Roman','Interpreter','latex','FontSize',7.5)


figure(51);
for kk=1:size(BB,1)
    plot(T,BB(kk,:)); hold on;
end

end