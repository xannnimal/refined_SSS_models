%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
Lin = 8;

%objFunc =  optimize_sensing_direction(angle); %objective function
% init = [-30 0]; %initial point
% lb = [-64 -64]; %lower bound
% ub = [64 64]; %upper  
x0 = zeros(144,1);
objFun = @(angle) optimize_sensing_direction(angle,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin);
lb = -1*ones(144,1);
ub = ones(144,1);
options = optimoptions('simulannealbnd','Display','iter');
[x,fval] = simulannealbnd(objFun,x0,lb,ub,options);
% zeros(144,3,'double'))