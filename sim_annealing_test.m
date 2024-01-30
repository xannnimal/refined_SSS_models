%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);

angle=phi_hat;
Lin = 8;
%objFunc =  optimize_sensing_direction(angle); %objective function
% init = [-30 0]; %initial point
% lb = [-64 -64]; %lower bound
% ub = [64 64]; %upper bound
opt = optimize_sensing_direction(angle,opm_matrix,R_hat,theta_hat,ch_types,Lin);
objFun = @opt;
x = simulannealbnd(objFun,angle);