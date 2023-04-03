f1=figure(1)
X=[1 2 3 4 5 6 7 8 9 10];
Y=[100 98 96 92 88 80 75 69 64.5 58; 72 47 35 22 15 10 7 4 3 2];
bar(X,Y); hold on;
% plot(X,Y,'-','LineWidth',1.5)
% title('Success Rate Vs Number of Obstacles in Randomized Enviroments','Interpreter','latex');
ylabel('Success Rate $(\%)$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
xlabel('Number of Obstacles','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
legend('DBaS-DDP','Penalty-DDP','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
% legend('DBaS-DDP','Penalty-DDP','FontName','Times New Roman');
f1.Position = [100  40  700  600];
%%
f2=figure(2)
x=[1 2 3 4 5 6 7 8 9 10];

y=100*[1 .996 .989 .985 .976 .954 .936 .917 .895 .843; .961 .924 .885 .825 .791 .75 .671 .651 .63 .566; .989 .973 .935 .885 .848 .782 .695 .65 .581 .524];

bar(x,y)
ylabel('Success Rate $(\%)$','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
xlabel('Number of Obstacles','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');
legend('DBaS-DDP','Penalty-DDP','CBF-DDP','FontName','Times New Roman','Fontsize',16,'Interpreter','latex');

% f2.Position = [699.6667  -42.3333  775.3333  660.0000];
f2.Position = [100  40  700  600];
