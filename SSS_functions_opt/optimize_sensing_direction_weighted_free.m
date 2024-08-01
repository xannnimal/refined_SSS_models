function conditionNumber = optimize_sensing_direction_weighted_free(angle_weight,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin)
% input: 144 angles
% fixed sensor positions calculate SSS matrix
% find condition number
% minimize condition number by changing angles

%Sin_vsh_vv: r_sphere,R,EX,EY,EZ,ch_types,Lin

angle = angle_weight(:,1);
weight = angle_weight(:,2);
phi_hat(:,2)  = weight .* angle;
[~,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',phi_hat',ch_types,Lin); %calculating SSS matrix

conditionNumber = cond(SNin);
end