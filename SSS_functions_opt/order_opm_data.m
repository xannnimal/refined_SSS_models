function unique_ordered_opm = order_opm_data(opm_matrix)

opm_matrix_x = opm_matrix(:,1);
opm_matrix_y = opm_matrix(:,2);
opm_matrix_z = opm_matrix(:,3);

count = 1;
num = 1;
dist_arr = zeros(144,1);
ordered_opm = zeros(144,3);

d = 0.0180;
d_hyp = sqrt(d^2 + d^2);
d_hyp = round(d_hyp, 3,'significant');


for j=1:length(opm_matrix)
    opm_mat_1 = opm_matrix(j,:);
    for i=1:length(opm_matrix)
        opm_mat_2 = opm_matrix(i,:);
        opm_mat_arr = [opm_mat_1;opm_mat_2];
        dist = pdist(opm_mat_arr);
        dist_arr(num) = dist;
        num = num + 1;
        dist_sig = round(dist, 3,'significant');
        if  dist_sig == d || dist_sig == d_hyp || dist_sig == 0
          ordered_opm(count,:) = opm_mat_2;
          count = count + 1;
        end
        
    end
end

unique_ordered_opm = unique(ordered_opm,'rows','stable');
% ordered_opm_matrix_x = unique_ordered_opm(:,1);
% ordered_opm_matrix_y = unique_ordered_opm(:,2);
% ordered_opm_matrix_z = unique_ordered_opm(:,3);
% 
% colormap = jet(length(opm_matrix));
% 
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
