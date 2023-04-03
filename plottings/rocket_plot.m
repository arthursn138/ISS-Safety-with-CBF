%% rocket plotting
function rocket_plot(T,X)
% T=0:dt:dt*(N-1);
% figure(1)
% subplot(1,2,1)
% plot(x0(2),x0(1),'or');
% hold on
% plot(xf(2),xf(1),'og');
% plot(X(2,:),X(1,:));
% figure(1)
% subplot(1,2,2);
% plot(T(1:end-1),U)
% subplot(1,3,3);
% plot(T,X(end,:))

% figure(2)
% plot3(x0(2),x0(3),x0(1),'or'); hold on;
% plot3(xf(2),xf(3),xf(1),'og');
% plot3(X(2,:),X(3,:),X(1,:));
% axis([-3 3 -3 3 0 4]); grid on
%%
% angle_gamma=[];
% angle_theta=[];
% Tnorm=[];
% angle_delta=[];
% w=[];
% for ii=1:length(X)-1
%     
%     w(ii)= norm([X(11,:);X(12,:);X(13,:)]);
% %     const2(ii) = X(1,ii) - tand(gamma_gs) * norm([X(2,ii);X(3,ii)]);
%     angle_gamma(ii) = atand(X(1,ii)/norm([X(2,ii);X(3,ii)]));
% 
% %     const3(ii) = 1 - 2 * (X(9,ii)^2 + X(10,ii)^2) - cosd(theta_max);
%     angle_theta(ii) = acosd(1 - 2 * (X(9,ii)^2 + X(10,ii)^2));
%     
%     Tnorm(ii) = norm([X(14,ii);X(15,ii);X(16,ii)]);
%     
%     angle_delta(ii) = acosd(X(14,ii)/norm([X(14,ii);X(15,ii);X(16,ii)]));
% end
% w(ii+1)=w(ii);
% angle_gamma(ii+1)=angle_gamma(ii);
% angle_theta(ii+1)=angle_theta(ii);
% Tnorm(ii+1)=Tnorm(ii);
% angle_delta(ii+1)=angle_delta(ii);
% 
% 
% figure(5)
% subplot(5,1,1)
% plot(T,w,T,w_max*ones(length(T))); ylabel('w');
% subplot(5,1,2)
% plot(T,angle_gamma,T,gamma_gs*ones(length(T))); ylabel('\gamma');
% subplot(5,1,3)
% plot(T,angle_theta,T,theta_max*ones(length(T))); ylabel('\theta');
% subplot(5,1,4)
% plot(T,Tnorm,T,Tmax*ones(length(T))); ylabel('T');
% subplot(5,1,5)
% plot(T,angle_delta,T,delta_max*ones(length(T))); ylabel('\delta');


% figure()
% plot(T,X(11,:), T, X(12,:), T, X(13,:))
end