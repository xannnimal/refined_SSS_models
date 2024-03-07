%% opm geometry from Peter at SANDIA
clear;
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
Lin = 1;
x0 = zeros(144,1);
sizeArr = 1;

% tic
% [Sin_test_old ,~]= Sin_vsh_vv([0,0,0]',opm_matrix(1:sizeArr,1:end)',R_hat(1:sizeArr,1:end)',theta_hat(1:sizeArr,1:end)',phi_hat(1:sizeArr,1:end)',ch_types(1:sizeArr),Lin);
% [Sin_test_new ,~]= Sin_vsh_vv_opt([0,0,0]',opm_matrix(1:sizeArr,1:end)',R_hat(1:sizeArr,1:end)',theta_hat(1:sizeArr,1:end)',phi_hat(1:sizeArr,1:end)',ch_types(1:sizeArr),Lin);
% toc

arr = zeros(1,8);
arr2 = zeros(1,8);
l_len = 2;
count = 1;
for l=1:l_len
    for m=-l:l
        arr(1,count) = prod(1:(l+m));
        count = count + 1;
    end
end

l_len = 2;
for l=1:l_len
    m=-l:l;
    arr2(1,l^2:length(m)+l^2 - 1) = prod(1:(l+m));
end

arr
arr2