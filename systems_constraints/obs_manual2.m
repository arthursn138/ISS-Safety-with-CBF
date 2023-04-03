function [h]=obs_manual2()
%     for ii=1:1
%         h{ii}= @(x) (x(1)^2+x(2)^2-1)^3-5*x(1)^2*x(2)^3;
%     end


% (((x)/0.5)^2+((y)/0.5)^2-1)^3 - 5*((x)/0.5)^2*((y)/0.5)^3


    h{1} = @(x) ((2*x(1))^2+(2*x(2))^2-1)^3-5*(2*x(1))^2*(2*x(2))^3; 
    h{2} = @(x) (( (((x(1)+2)*cos(90)-(x(2)+2)*sin(90))/0.5)^2+(((x(2)+2)*cos(90)+(x(1)+2)*sin(90))/0.5)^2-1)^3 - 5*(((x(1)+2)*cos(90)-(x(2)+2)*sin(90))/0.5)^2*(((x(2)+2)*cos(90)+(x(1)+2)*sin(90))/0.5)^3); 
    h{3} = @(x) (5*(x(1)-0.5)^2 + (x(2)+3)^2) -1/2; 
    h{4} = @(x) abs(x(1)+2)+abs(x(2)-2) - 1; 
    h{5} = @(x) abs((x(1)-3)+(x(2)+1))+abs((x(1)-3)-(x(2)+1)) - 1; 
    h{6} = @(x) (5*(x(1)-1)^2 + (x(2)+2)^2) -1; 
    h{7} = @(x) ((x(1)+0.5)^2 + 50*(x(2)+4)^2) -2; 
    h{8} = @(x) ((x(1)-1.5)^2 + (x(2)-2.5)^2) -1; 
    
    figure(1) %heart 1
    h1=fimplicit(@(x,y) ((2*x)^2+(2*y)^2-1)^3-5*(2*x)^2*(2*y)^3 ,[-2,2,-2,2]);
    xd1=h1.XData; yd1=h1.YData;
%     fill(xd1,yd1,'r'); 
    fill(xd1,yd1,[0.6350 0.0780 0.1840]);
%     patch([xd1],[yd1],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd1],[yd1],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on
    
    %heart 2
    h2=fimplicit(@(x,y) (( (((x+2)*cos(90)-(y+2)*sin(90))/0.5)^2+(((y+2)*cos(90)+(x+2)*sin(90))/0.5)^2-1)^3 - 5*(((x+2)*cos(90)-(y+2)*sin(90))/0.5)^2*(((y+2)*cos(90)+(x+2)*sin(90))/0.5)^3),[-3,0,-3,-1]);
    xd2=h2.XData; yd2=h2.YData;
    fill(xd2,yd2,[0.6350 0.0780 0.1840]);
%     patch([xd2],[yd2],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd2],[yd2],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on

    % circle
    h3=fimplicit(@(x,y) (5*(x-.5)^2 + (y+3)^2) -1/2 ,[-1,3,-4.1,-.9]);
    xd3=h3.XData; yd3=h3.YData;
    fill(xd3,yd3,[0.6350 0.0780 0.1840]);
%     patch([xd3],[yd3],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd3],[yd3],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on
    
    % circle
    h4=fimplicit(@(x,y) (5*(x-1)^2 + (y+2)^2) -1 ,[-1,3,-3.1,-.9]);
    xd4=h4.XData; yd4=h4.YData;
    fill(xd4,yd4,[0.6350 0.0780 0.1840]);
%     patch([xd4],[yd4],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    patch([xd4],[yd4],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on
    
    % circle
    h5=fimplicit(@(x,y) ((x+0.5)^2 + 50*(y+4)^2) -2 ,[-5,3,-6.1,0]);
    xd5=h5.XData; yd5=h5.YData;
%     patch([xd5],[yd5],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    fill(xd5,yd5,[0.6350 0.0780 0.1840]);
%     patch([xd5],[yd5],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on
    
    % circle
    h6=fimplicit(@(x,y) ((x-1.5)^2 + (y-2.5)^2) -1 ,[0,5,0,4]);
    xd6=h6.XData; yd6=h6.YData;
%     patch([xd6],[yd6],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);    
    fill(xd6,yd6,[0.6350 0.0780 0.1840]);
%     patch([xd6],[yd6],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on    
    
    % diamond
    h7=fimplicit(@(x,y) abs(x+2)+abs(y-2) - 1,[-3,1,1,3])  ;  
    xd7=h7.XData; yd7=h7.YData;
%     patch([xd7],[yd7],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    fill(xd7,yd7,[0.6350 0.0780 0.1840]);
%     patch([xd7],[yd7],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on;
    
    % square
    h8=fimplicit(@(x,y) abs((x-3)+(y+1))+abs((x-3)-(y+1)) - 1,[2,4,-1.6,0])  ;  
    xd8=h8.XData; yd8=h8.YData;
%     patch([xd8],[yd8],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);    
    fill(xd8,yd8,[0.6350 0.0780 0.1840]);
    patch([xd8],[yd8],[0.6350 0.0780 0.1840],'FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
    hold on;
    
    box on;
    
%     figure(1000)
%     patch([xd1],[yd1],'r','FaceAlpha',
%     1,'linestyle','none','EdgeAlpha',0); hold on
%     patch([xd2],[yd2],'y','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd3],[yd3],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd4],[yd4],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd5],[yd5],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd6],[yd6],'r','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd7],[yd7],'b','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);
%     patch([xd8],[yd8],'b','FaceAlpha', 1,'linestyle','none','EdgeAlpha',0);

%     axis([-4.5,4.5,-4.5,4.5]);
%     axis equal


%     h4=fimplicit(@(x,y) 1/(1-(x-2)^2) + 1/(1-y^2) - 5,[1,3,-.9,.9]);
%     h4=fimplicit(@(x,y) ( 1/(1-x^2) + 1/(1-y^2)) - 4,[-.9,.9,-.9,.9])  ;  
%     xd=h4.XData; yd=h4.YData;
%     fill(xd,yd,'b');
%     hold on;

%    obs_plot=[circ1;circ2;circ3;circ4];
%    figure(10000)
%    for ii=1:num_obs
%    rectangle('Position',circ(ii,:),'Curvature',[1 1],'FaceColor',[0.3010 0.7450 0.9330]); hold on
%    end
% clear x y
% x=-1.50274:0.001:1.50274;
% 
%    for ii=1:length(x)
%     sol = roots([1, 0, (3*x(ii)^2 - 3), - 5*x(ii)^2, (- 6*x(ii)^2 + 3*x(ii)^4 + 3), 0, x(ii)^6 - 3*x(ii)^4 + 3*x(ii)^2 - 1]);
%     sol=sol(imag(sol)==0);
%     y(1,ii)=sol(1);
%     sol = sol(sol~=sol(1));
%     y(2,ii)=sol(1);
%    end
% 
%    figure(1000)
%    plot(x,y); hold on;
%    fill(x,y,'r')

% ezplot('(x^2+y^2-1)^3 = x^2*y^3',[-1.5,1.5,-1.2,1.8])
clc

end
    
