%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry_avg(filename);
Lin = 8;
sensor_len = length(opm_matrix);
x0_arr = zeros(sensor_len,2);

objFun = @(angle_weight) optimize_sensing_direction_weighted(angle_weight,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin);
lb_weight = [-pi*ones(sensor_len,1),0.1*ones(sensor_len,1)]; %restrict to positive numbers, 0.1 to 10
ub_weight = [pi*ones(sensor_len,1),10*ones(sensor_len,1)];
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',1000,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
[angle_weight,fval] = simulannealbnd(objFun,x0_arr,lb_weight,ub_weight,options);
angles = angle_weight(:,1);
weight = angle_weight(:,2);

sensing_dir = cos(angles).*phi_hat + sin(angles).*theta_hat;
