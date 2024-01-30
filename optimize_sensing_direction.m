function conditionNumber = optimize_sensing_direction(angle,opm_matrix,R_hat,theta_hat,ch_types,Lin)
% input: 144 angles
% fixed sensor positions calculate SSS matrix
% find condition number
% minimize condition number by changing angles
% where should I put the angle? Theta or phi?

%Sin_vsh_vv: r_sphere,R,EX,EY,EZ,ch_types,Lin
%opt phi_hat
[~,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',angle',ch_types,Lin); %calculating SSS matrix
conditionNumber = cond(SNin);
end