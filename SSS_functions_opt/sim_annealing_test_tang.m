%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
% [opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry_avg(filename);
Lin = 8;
sensor_len = length(opm_matrix);
x0 = phi_hat(:,2);

objFun_tang = @(angles_tang) optimize_sensing_direction_tang(angles_tang,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin);
lb = -1*ones(sensor_len,1);
ub = ones(sensor_len,1);
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',1000,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
[angles_tang,fval_tang] = simulannealbnd(objFun_tang,x0,lb,ub,options);

phi_hat_new = [phi_hat(:,1), angles_tang, phi_hat(:,3)];

% subspace_angles = zeros(144,1);

% for i=(1:length(phi_hat))
%     subspace_angles(i)=subspace(phi_hat_new(i,:)',phi_hat(i,:)');
% end
% 
% histogram(subspace_angles*180/pi,5)
% title('Deviation between initial and optimized sensing direction')
% xlabel('Subspace angle (degrees)')
% ylabel('Count')
% 
% figure(2);
% hold on
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'DisplayName','Data')
% q1 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), phi_hat(:,1), phi_hat(:,2), phi_hat(:,3),'DisplayName','Original Sensing Vector');
% q1.Color = "#0072BD";
% q2 = quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3), phi_hat_new(:,1), phi_hat_new(:,2), phi_hat_new(:,3),'DisplayName','Optimized Sensing Vector (unrestricted)');
% q2.Color = "#7E2F8E";
% grid on
% rotate3d
% view(135, 20);
% hold off
