clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
len=144;
% condition_num_arr= zeros(len,1);
condition_num_arr = load("condition_num_var_arr.mat");
condition_num_arr = condition_num_arr.condition_num_arr;
Lin = 8;

% for j = 1:length(condition_num_arr)
%    [opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry_var(filename,j);
%    [~,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin);
%     condition_num_arr(j) = cond(SNin);
% end
xarr = linspace(1,len,len);
scatter(xarr,condition_num_arr,"filled")
% ylim([0,1.4*1e5])
