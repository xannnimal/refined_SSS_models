%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
unique_ordered_opm = order_opm_data(opm_matrix);

[opm_matrix_avg,R_hat_avg,theta_hat_avg2,phi_hat_avg2,~] = gen_opm_geometry_avg_full(filename);

Lin = 8;
sensor_len = length(opm_matrix);
x0 = zeros(sensor_len,1);

phi_hat_avg = zeros(sensor_len,3);
theta_hat_avg = zeros(sensor_len,3);
% 
% for j=1:(sensor_len/4)
%     if j==1
%         avg_phi = mean(phi_hat(j:4*j,:));
%         avg_theta = mean(theta_hat(j:4*j,:));
%         phi_hat_avg(j:4*j,:) = repmat(avg_phi,4,1);
%         theta_hat_avg(j:4*j,:) = repmat(avg_theta,4,1);
%     else
%         avg_phi = mean(phi_hat(4*j-3:4*j,:));
%         avg_theta = mean(theta_hat(4*j-3:4*j,:));
%         phi_hat_avg(4*j-3:4*j,:) = repmat(avg_phi,4,1);
%         theta_hat_avg(4*j-3:4*j,:) = repmat(avg_theta,4,1);
%     end
% end

objFun = @(angles) optimize_sensing_direction(angles,unique_ordered_opm,R_hat,theta_hat_avg2,phi_hat_avg2,ch_types,Lin);
lb = -pi*ones(sensor_len,1);
ub = pi*ones(sensor_len,1);
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',50, ...
    'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplottemperature},'AnnealingFcn', @annealingfast_avg);
[angles,fval,temperature] = simulannealbnd(objFun,x0,lb,ub,options);
sensing_dir = cos(angles).*phi_hat_avg2+ sin(angles).*theta_hat_avg2;



