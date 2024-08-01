filename="headwithsensors1.mat";
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);

coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";
[R,EX,EY,EZ, ch_types] = gen_squid_geometry(rawfile,coordsys);
EZ = EZ';
xdir = EZ(:,1);
ydir = EZ(:,2);
zdir = EZ(:,3);
r = sqrt(sum(EZ.^2,2)); %already normalized, just precaution

theta=r .* acos(zdir./r);
phi= sign(ydir).*acos(xdir./sqrt(xdir.^2+ ydir.^2));

for i=(1:size(xdir,1))
    theta2(i)=acos(zdir(i)/sqrt(xdir(i)^2+ ydir(i)^2+ zdir(i)^2));
    phi2(i)= sign(ydir(i))*acos(xdir(i)/sqrt(xdir(i)^2+ ydir(i)^2));

end
theta2 = theta2';
phi2 = phi2';
% nvar=36;
% y = randn(nvar,1);
% y = y./norm(y);
% y_144 = nvar*4;
% 
% for j=1:nvar
%     if j==1
%         y_144(j:4*j,:) = repmat(y(j),4,1);
%     else
%         y_144(4*j-3:4*j,:) = repmat(y(j),4,1);
%     end
% end
% newx = y_144;

% opm_matrix_x = opm_matrix(:,1);
% opm_matrix_y = opm_matrix(:,2);
% opm_matrix_z = opm_matrix(:,3);
% 
% count = 1;
% num = 1;
% dist_arr = zeros(144,1);
% ordered_opm = zeros(144,3);
% near = pdist(opm_matrix(1:2,:));
% 
% d = 0.0180;
% d_hyp = sqrt(d^2 + d^2);
% d_hyp = round(d_hyp, 3,'significant');
% 
% 
% for j=1:length(opm_matrix)
%     opm_mat_1 = opm_matrix(j,:);
%     for i=1:length(opm_matrix)
%         opm_mat_2 = opm_matrix(i,:);
%         opm_mat_arr = [opm_mat_1;opm_mat_2];
%         dist = pdist(opm_mat_arr);
%         dist_arr(num) = dist;
%         num = num + 1;
%         dist_sig = round(dist, 3,'significant');
%         if  dist_sig == d || dist_sig == d_hyp || dist_sig == 0
%           ordered_opm(count,:) = opm_mat_2;
%           count = count + 1;
%         end
% 
%     end
% end
% 
% unique_ordered_opm = unique(ordered_opm,'rows','stable');
% ordered_opm_matrix_x = unique_ordered_opm(:,1);
% ordered_opm_matrix_y = unique_ordered_opm(:,2);
% ordered_opm_matrix_z = unique_ordered_opm(:,3);
% 
% colormap = jet(length(opm_matrix));

  
% figure(4);
% hold on
% for k = 1:6
%     if k==1
%         scatter3(opm_matrix_x(1:4),opm_matrix_y(1:4*k),opm_matrix_z(1:4*k),'filled','Color',colormap(k,:));
%     else
%         scatter3(opm_matrix_x(4*k-3:4*k),opm_matrix_y(4*k-3:4*k),opm_matrix_z(4*k-3:4*k),'filled','Color',colormap(k,:));
%     end
% end
% grid on
% rotate3d
% view(135, 20);
% hold off
% 
% figure(5);
% hold on
% for k = 1:6
%     if k==1
%         scatter3(ordered_opm_matrix_x(1:4),ordered_opm_matrix_y(1:4*k),ordered_opm_matrix_z(1:4*k),'filled','Color',colormap(k,:));
%     else
%         scatter3(ordered_opm_matrix_x(4*k-3:4*k),ordered_opm_matrix_y(4*k-3:4*k),ordered_opm_matrix_z(4*k-3:4*k),'filled','Color',colormap(k,:));
%     end
% end
% grid on
% rotate3d
% view(135, 20);
% hold off

% 
% figure(6);
% hold on
% scatter3(opm_matrix_x(1:4),opm_matrix_y(1:4),opm_matrix_z(1:4),'red','filled');
% scatter3(opm_matrix_x(5:8),opm_matrix_y(5:8),opm_matrix_z(5:8),'blue','filled');
% scatter3(opm_matrix_x(9:12),opm_matrix_y(9:12),opm_matrix_z(9:12),'green','filled');
% scatter3(opm_matrix_x(13:16),opm_matrix_y(13:16),opm_matrix_z(13:16),'magenta','filled');
% scatter3(opm_matrix_x(17:20),opm_matrix_y(17:20),opm_matrix_z(17:20),'cyan','filled');
% scatter3(opm_matrix_x(21:24),opm_matrix_y(21:24),opm_matrix_z(21:24),'filled','markerfacecolor','#EDB120');
% scatter3(opm_matrix_x(25:28),opm_matrix_y(25:28),opm_matrix_z(25:28),'filled','markerfacecolor','#77AC30');
% scatter3(opm_matrix_x(29:32),opm_matrix_y(29:32),opm_matrix_z(29:32),'filled','markerfacecolor','#7E2F8E');
% grid on
% rotate3d
% view(135, 20);
% hold off

% figure(7);
% hold on
% scatter3(opm_matrix_x(1),opm_matrix_y(1),opm_matrix_z(1),'red','filled','DisplayName','p1');
% scatter3(opm_matrix_x(2),opm_matrix_y(2),opm_matrix_z(2),'blue','filled','DisplayName','p2');
% scatter3(opm_matrix_x(3),opm_matrix_y(3),opm_matrix_z(3),'green','filled','DisplayName','p3');
% scatter3(opm_matrix_x(4),opm_matrix_y(4),opm_matrix_z(4),'magenta','filled','DisplayName','p4');
% grid on
% rotate3d
% view(135, 20);
% hold off


% %% opm geometry from Peter at SANDIA
% clear;
% 
% coordsys = 'device'; 
% rawfile = "sample_audvis_raw.fif";
% 
% [R,EX,EY,EZ, ch_types] = gen_squid_geometry(rawfile,coordsys);
% R = R';
% EX = EX';
% EY = EY';
% EZ = EZ';
% Lin = 8;
% 
% sensor_len = length(R);
% x0 = zeros(sensor_len,1);
% x0_arr = zeros(sensor_len,2);

% filename="headwithsensors1.mat";
% %generate helmet pos and ori with "gen_opm_geometry"
% [opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
% Lin = 1;
% x0 = zeros(144,1);
% sizeArr = 1;
% 
% weights = load('weighted_weights_1000.mat');
% weights = weights.weight;
% 
% histogram(weights)

% weights_only = load("squid_weights_only_1000.mat");
% weights_only = weights_only.weights;
% sensing_dir_weighted = weights_only .* EY;
% 
% objFun = @(weight) optimize_sensing_direction_weighted_only(weight,R,EX,sensing_dir_weighted,EZ,ch_types,Lin);
% lb_weight = 0.1*ones(sensor_len,1); %restrict to positive numbers, 0.1 to 10
% ub_weight = 10*ones(sensor_len,1);
% options = optimoptions('simulannealbnd','Display','iter','MaxIterations',1000,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
% [weights,fval] = simulannealbnd(objFun,x0,lb_weight,ub_weight,options);

% tic
% [Sin_test_old ,~]= Sin_vsh_vv([0,0,0]',opm_matrix(1:sizeArr,1:end)',R_hat(1:sizeArr,1:end)',theta_hat(1:sizeArr,1:end)',phi_hat(1:sizeArr,1:end)',ch_types(1:sizeArr),Lin);
% [Sin_test_new ,~]= Sin_vsh_vv_opt([0,0,0]',opm_matrix(1:sizeArr,1:end)',R_hat(1:sizeArr,1:end)',theta_hat(1:sizeArr,1:end)',phi_hat(1:sizeArr,1:end)',ch_types(1:sizeArr),Lin);
% toc

% arr = zeros(1,8);
% arr2 = zeros(1,8);
% l_len = 1;
% count = 1;
% for l=1:l_len
%     for m=-l:l
%         arr(1,count) = prod(1:(l+m));
%         count = count + 1;
%     end
% end
% 
% l_len = 2;
% for l=1:l_len
%     m=-l:l;
%     arr2(1,l^2:length(m)+l^2 - 1) = prod(1:(l+m));
% end
% 
% arr;
% arr2

% reshape(bsxfun(@plus,A,(-3:2).'),1,[])
