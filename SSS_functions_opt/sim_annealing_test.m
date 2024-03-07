%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
Lin = 8;
x0 = zeros(144,1);

objFun = @(angle) optimize_sensing_direction(angle,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin);
lb = -1*ones(144,1);
ub = ones(144,1);
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',1000,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
[phi_hat_y,fval] = simulannealbnd(objFun,x0,lb,ub,options);
% phi_hat_new = [phi_hat(:,1), phi_hat_y, phi_hat(:,3)];

phi_hat_data = load('phi_hat_opt.mat');
phi_hat_new = phi_hat_data.phi_hat_opt;

subspace_angles = zeros(144,1);

for i=(1:length(phi_hat))
    subspace_angles(i)=subspace(phi_hat_new(i,:)',phi_hat(i,:)');
end

histogram(subspace_angles*180/pi,5)
title('Deviation between initial and optimized sensing direction')
xlabel('Subspace angle (degrees)')
ylabel('Count')

figure(7);
hold on
scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'DisplayName','Data')
q1 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), phi_hat(:,1), phi_hat(:,2), phi_hat(:,3),'DisplayName','Original Sensing Vector');
q1.Color = "#0072BD";
q2 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), phi_hat_new(:,1), phi_hat_new(:,2), phi_hat_new(:,3),'DisplayName','Optimized Sensing Vector');
q2.Color = "#D95319";
grid on
rotate3d
view(135, 20);
hold off