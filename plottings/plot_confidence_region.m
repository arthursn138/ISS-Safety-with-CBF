% plotting confidence region
function plot_confidence_region(fig_num,X1,X2,X3,N,zstar,spc,facecolor,alpha,edgecolor,linestyle)
    x1m=mean(X1); x2m=mean(X2); x3m=mean(X3); num_samples=size(X1,1);
    for ii=1:spc:N
        x=X1(:,ii); y=X2(:,ii); z=X3(:,ii);
        CovXY=(x-mean(x))'*(y-mean(y))/(num_samples-1);
        CovYX=CovXY;    
        CovXZ=(x-mean(x))'*(z-mean(z))/(num_samples-1);
        CovZX=CovXZ;
        CovYZ=(y-mean(y))'*(z-mean(z))/(num_samples-1);
        CovZY=CovYZ;
        covmatrix3d=[var(x) CovXY CovXZ;CovYX var(y) CovYZ;CovZX CovZY var(z)];
       
        cov_matrix3d_regularized=covmatrix3d*zstar;
        [a,b]=eig(cov_matrix3d_regularized,'vector');
        [b,ind] = sort(b,'descend');
        a = a(:, ind);
        v=[sqrt(b(1))*a(:,1) sqrt(b(2))*a(:,2) sqrt(b(3))*a(:,3)];
        mean3D = [x1m(ii); x2m(ii); x3m(ii)];

        [x,y,z] = sphere(round(sqrt(N-1)));

        surfmatrix=[x(:),y(:),z(:)]';
        extended3dellip = real(v)*surfmatrix + repmat(mean3D, 1, size(surfmatrix,2));
 
        x = reshape(extended3dellip(1,:),size(x));
        y = reshape(extended3dellip(2,:),size(y));
        z = reshape(extended3dellip(3,:),size(z));
        figure(fig_num);
        conf_reg_3d = surf(x,y,z);
        conf_reg_3d.EdgeColor = edgecolor;
        conf_reg_3d.LineStyle = linestyle;
        conf_reg_3d.FaceColor = facecolor;
        conf_reg_3d.FaceAlpha = alpha;
        hold on

    
%     x = v(1,:)*[cos(t).*cos(u);sin(t).*cos(u);sin(u)] + x1m(ii);
%     y = v(2,:)*[cos(t).*cos(u);sin(t).*cos(u);sin(u)] + x2m(ii);
%     z = v(3,:)*[cos(t).*cos(u);sin(t).*cos(u);sin(u)] + x3m(ii);
%     fill3([x],[y],[z],color,'FaceColor','none','FaceAlpha', alpha,'EdgeAlpha',1);

% to check with some exsisting work in matlab exchange
% https://www.mathworks.com/matlabcentral/fileexchange/16543-plot_gaussian_ellipsoid
    end
end