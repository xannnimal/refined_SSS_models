clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
[~,~,~,phi_hat_avg2,~] = gen_opm_geometry_avg_full(filename);


sensing_dir_opt = load("grouped_sensing_direction_5000.mat");
sensing_dir_opt = sensing_dir_opt.sensing_dir;
ordered_sensing_dir_opt = load("ordered_dir_1000.mat");

ordered_dir_opt = ordered_sensing_dir_opt.sensing_dir;

ordered_dir_opt_200 = load("ordered_dir_200.mat");
ordered_dir_opt_200 = ordered_dir_opt_200.sensing_dir;

opt_dir_1000= load("sensing_dir_1000.mat");
opt_dir_1000 = opt_dir_1000.sensing_dir;

ordered_opm = order_opm_data(opm_matrix);

opm_x = ordered_opm(:,1);
opm_y = ordered_opm(:,2);
opm_z = ordered_opm(:,3);
og_dir_x = phi_hat_avg2(:,1);
og_dir_y = phi_hat_avg2(:,2);
og_dir_z = phi_hat_avg2(:,3);
opt_dir_x = opt_dir_1000(:,1);
opt_dir_y = opt_dir_1000(:,2);
opt_dir_z = opt_dir_1000(:,3);

n=1;
m=144;
count_arr=zeros(length(opm_matrix),1);

% for j=1:length(opm_matrix)
%     dir1 = ordered_dir_opt(j,:);
%     count=0;
%     for i=1:length(opm_matrix)
%         dir2 = ordered_dir_opt(i,:);
%         if isequal(dir1, dir2)
%           count =  count + 1;
%         end
%     end
%     count_arr(j) = count;
% end

% figure(6);
% hold on
% scatter3(opm_x(n:m),opm_y(n:m),opm_z(n:m),'DisplayName','Data');
% q1 = quiver3(opm_x(n:m),opm_y(n:m),opm_z(n:m), og_dir_x(n:m),og_dir_y(n:m),og_dir_z(n:m),'DisplayName','Original Sensing Vector');
% q1.Color = "#0072BD";
% % scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'DisplayName','Data')
% q2 = quiver3(opm_x(n:m),opm_y(n:m),opm_z(n:m),opt_dir_x(n:m),opt_dir_y(n:m),opt_dir_z(n:m),'DisplayName','Optimized Sensing Vector');
% q2.Color = "#D95319";
% grid on
% rotate3d
% view(135, 20);
% hold off

figure(7);
hold on
scatter3(opm_x,opm_y,opm_z,'DisplayName','Data');
q1 = quiver3(opm_x,opm_y,opm_z, og_dir_x,og_dir_y,og_dir_z,'DisplayName','Phi hat');
q1.Color = "#0072BD";
q2 = quiver3(opm_x,opm_y,opm_z,theta_hat(:,1),theta_hat(:,2),theta_hat(:,3),'DisplayName','Theta hat');
q2.Color = "#D95319";
q3 = quiver3(opm_x,opm_y,opm_z,opt_dir_x,opt_dir_y,opt_dir_z,'DisplayName','Optimized Sensing Vector');
q3.Color = "#d919bc";
grid on
rotate3d
view(135, 20);
hold off


% figure(4);
% hold on
% scatter3(opm_matrix_x(1:4),opm_matrix_y(1:4),opm_matrix_z(1:4),'red','filled');
% scatter3(opm_matrix_x(4:8),opm_matrix_y(4:8),opm_matrix_z(4:8),'blue','filled');
% scatter3(opm_matrix_x(8:12),opm_matrix_y(8:12),opm_matrix_z(8:12),'green','filled');
% scatter3(opm_matrix_x(12:16),opm_matrix_y(12:16),opm_matrix_z(12:16),'magenta','filled');
% scatter3(opm_matrix_x(16:20),opm_matrix_y(16:20),opm_matrix_z(16:20),'cyan','filled');
% scatter3(opm_matrix_x(20:24),opm_matrix_y(20:24),opm_matrix_z(20:24),'filled','markerfacecolor','#EDB120');
% scatter3(opm_matrix_x(24:28),opm_matrix_y(24:28),opm_matrix_z(24:28),'filled','markerfacecolor','#77AC30');
% scatter3(opm_matrix_x(28:32),opm_matrix_y(28:32),opm_matrix_z(28:32),'filled','markerfacecolor','#7E2F8E');
% grid on
% rotate3d
% view(135, 20);
% hold off


