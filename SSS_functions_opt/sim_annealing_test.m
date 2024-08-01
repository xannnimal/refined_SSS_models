%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry_avg(filename);
[opm_matrix_full,R_hat_full,theta_hat_full,phi_hat_full,ch_types_full] = gen_opm_geometry(filename);
[opm_matrix_avg_144,R_hat_avg_144,theta_hat_avg_144,phi_hat_avg_144,ch_types_avg_144] = gen_opm_geometry_avg_144(filename);

Lin = 8;
sensor_len = length(opm_matrix);
sensor_len_full = length(opm_matrix_full);
x0 = zeros(sensor_len,1);
x0_full = zeros(sensor_len_full,1);
type = 1;
type_144 =  0;

% objFun = @(angles) optimize_sensing_direction_avg(angles,opm_matrix_avg_144,R_hat_avg_144,theta_hat_avg_144, ...
%     phi_hat_avg_144,ch_types_avg_144,Lin,type_144);
% lb = -pi*ones(sensor_len_full,1);
% ub = pi*ones(sensor_len_full,1);
% options = optimoptions('simulannealbnd','Display','iter','MaxIterations',100,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
% [angles_144,fval] = simulannealbnd(objFun,x0_full,lb,ub,options);

objFun_avg = @(angles_full) optimize_sensing_direction_avg(angles_full,opm_matrix_full,R_hat_full,phi_hat_full,theta_hat_full, ...
    ch_types_full,Lin,type);
lb_full = -pi*ones(sensor_len_full,1);
ub_full = pi*ones(sensor_len_full,1);
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',100,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
[angles_avg,fval_avg] = simulannealbnd(objFun_avg,x0_full,lb_full,ub_full,options);

% objFun = @(angles) optimize_sensing_direction(angles,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin,type);
% lb = -pi*ones(sensor_len_full,1);
% ub = pi*ones(sensor_len_full,1);
% options = optimoptions('simulannealbnd','Display','iter','MaxIterations',10,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
% [angles,fval] = simulannealbnd(objFun,x0,lb,ub,options);

% objFun_full = @(angles_full) optimize_sensing_direction(angles_full,opm_matrix_full,R_hat_full,phi_hat_full,theta_hat_full,ch_types_full,Lin);
% lb_full = -pi*ones(sensor_len,1);
% ub_full = pi*ones(sensor_len,1);
% options = optimoptions('simulannealbnd','Display','iter','MaxIterations',10,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
% [angles_full,fval_full] = simulannealbnd(objFun_full,x0_full,lb_full,ub_full,options);

% sensing_dir = cos(angles).*phi_hat + sin(angles).*theta_hat;

% sensing_dir = cos(angles).*phi_hat + sin(angles).*theta_hat;

% figure(7);
% hold on
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'DisplayName','Data')
% scatter3(sensing_dir(:,1),sensing_dir(:,2),sensing_dir(:,3),'DisplayName','Optimized')
% grid on
% rotate3d
% view(135, 20);
% hold off
% 
% figure(1);
% hold on
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'DisplayName','Data')
% q1 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), phi_hat(:,1), phi_hat(:,2), phi_hat(:,3),'DisplayName','Original Sensing Vector');
% q1.Color = "#0072BD";
% scatter3(opm_matrix_full(:,1),opm_matrix_full(:,2),opm_matrix_full(:,3),'DisplayName','Data')
% q2 = quiver3(opm_matrix_full(:,1),opm_matrix_full(:,2),opm_matrix_full(:,3), phi_hat_full(:,1), phi_hat_full(:,2), phi_hat_full(:,3),'DisplayName','Original Sensing Vector');
% q2.Color = "#D95319";
% 
% q2 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), sensing_dir(:,1), sensing_dir(:,2), sensing_dir(:,3),'DisplayName','Optimized Sensing Vector');
% q2.Color = "#D95319";
% grid on
% rotate3d
% view(135, 20);
% hold off
% 
% figure(2);
% hold on
% scatter3(opm_matrix_full(:,1),opm_matrix_full(:,2),opm_matrix_full(:,3),'DisplayName','Data')
% q1 = quiver3(opm_matrix_full(:,1),opm_matrix_full(:,2),opm_matrix_full(:,3), phi_hat_full(:,1), phi_hat_full(:,2), phi_hat_full(:,3),'DisplayName','Original Sensing Vector');
% q1.Color = "#D95319";
% q2 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), sensing_dir(:,1), sensing_dir(:,2), sensing_dir(:,3),'DisplayName','Optimized Sensing Vector');
% q2.Color = "#D95319";
% grid on
% rotate3d
% view(135, 20);
% hold off

