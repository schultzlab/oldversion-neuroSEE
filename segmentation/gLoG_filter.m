function [masks] = gLoG_filter(red_initial_image)

theta=0:45:135;
sigma_x = 4:8;
sigma_y = 4:8;
a = [];
b = [];
c = [];
lap_x = [];
lap_y = [];
G_x_y = [];
glog_kernel = [];
norm_glog_kernel = [];
lognorm_glog_kernel = [];
lambda = 1; %normalising factor

%Compute the parameters that control the shape and orientation of the gLoG
%kernels: a,b,c
for iii=1:4
    for ii=1:5
        for i=1:5
            a(i,ii,iii) = (cos(theta(iii))^2 /(2*sigma_x(i)^2)) + (sin(theta(iii))^2/(2*sigma_y(ii)^2));
            b(i,ii,iii) = (-sin(2*theta(iii))/(4*sigma_x(i)^2)) + (sin(2*theta(iii))/(4*sigma_y(ii)^2));
            c(i,ii,iii) = (sin(theta(iii))^2 /(2*sigma_x(i)^2)) + (cos(theta(iii))^2/(2*sigma_y(ii)^2));
        end
    end
end

% Reshape the arrays
a_ = reshape(a, [1 100]);
b_ = reshape(b, [1 100]);
c_ = reshape(c, [1 100]);

% Compute for G(x,y), the Laplacians and the gLoG kernels
for i=1:100
    for x=1:size(red_initial_image,1)
        for y=1:size(red_initial_image,2)
            G_x_y(x,y,i) = lambda*exp(-(a_(i)*x^2 + 2*b_(i)*x*y + c_(i)*y^2));
            lap_x(x,y,i) = ((2*a_(i)*x + 2*b_(i)*y)^2 - 2*a_(i))*G_x_y(x,y,i);
            lap_y(x,y,i) = ((2*b_(i)*x + 2*c_(i)*y)^2 - 2*c_(i))*G_x_y(x,y,i);
            glog_kernel(x,y,i) = lap_x(x,y,i) + lap_y(x,y,i);
        end
    end
end

% Scale-normalise the kernels 
counter = 1;
summed_gk_theta_norm = zeros(size(red_initial_image,1),size(red_initial_image,2),length(theta));
summed_gk_theta_lognorm = zeros(size(red_initial_image,1),size(red_initial_image,2),length(theta));
L_xy_norm = [];
L_xy_lognorm = [];
for th=1:4
    for sy=1:5
        for sx = 1:5
            norm_glog_kernel(:,:,counter) = sigma_x(sx)*sigma_y(sy)*glog_kernel(:,:,counter);
            lognorm_glog_kernel(:,:,counter) = (1+0.5*log(sigma_x(sx)))*(1+0.5*log(sigma_y(sy)))*glog_kernel(:,:,counter);
            counter = counter+1;
        end
    end
end
% Create the summed gLoG filter
summed_gk_theta_norm(:,:,1) = (1/25)*sum(norm_glog_kernel(:,:,1:25),3);
summed_gk_theta_norm(:,:,2) = (1/25)*sum(norm_glog_kernel(:,:,26:50),3);
summed_gk_theta_norm(:,:,3) = (1/25)*sum(norm_glog_kernel(:,:,51:75),3);
summed_gk_theta_norm(:,:,4) = (1/25)*sum(norm_glog_kernel(:,:,76:100),3);
summed_gk_theta_lognorm(:,:,1) = (1/25)*sum(lognorm_glog_kernel(:,:,1:25),3);
summed_gk_theta_lognorm(:,:,2) = (1/25)*sum(lognorm_glog_kernel(:,:,26:50),3);
summed_gk_theta_lognorm(:,:,3) = (1/25)*sum(lognorm_glog_kernel(:,:,51:75),3);
summed_gk_theta_lognorm(:,:,4) = (1/25)*sum(lognorm_glog_kernel(:,:,76:100),3);

% Construct the scale-space response maps by convolving the summed gLoG filter with the
% image 
for i=1:4
    L_xy_norm(:,:,i) = conv2(red_initial_image,summed_gk_theta_norm(:,:,i));
    L_xy_lognorm(:,:,i) = conv2(red_initial_image,summed_gk_theta_lognorm(:,:,i));
end
L_xy_norm = L_xy_norm(1:size(red_initial_image,1),1:size(red_initial_image,2),:);
L_xy_lognorm = L_xy_lognorm(1:size(red_initial_image,1),1:size(red_initial_image,2),:);



end