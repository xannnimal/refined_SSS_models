function conditionNumber = optimize_sensing_direction_avg(angle,opm_matrix,R_hat,phi_hat,theta_hat,ch_types,Lin,type)

if type==1
    phi_hat_avg_len = length(phi_hat);
    phi_hat_avg = zeros(phi_hat_avg_len,3);
    
    for j=1:(phi_hat_avg_len/4)
        if j==1
            avg = mean(phi_hat(j:4*j,:));
            phi_hat_avg(j:4*j,:) = repmat(avg,4,1);
        else
            avg = mean(phi_hat(4*j-3:4*j,:));
            phi_hat_avg(4*j-3:4*j,:) = repmat(avg,4,1);
        end
    end
    
    sensing_dir = cos(angle).*phi_hat_avg + sin(angle).*theta_hat;
    % phi_hat(:,2)= angle;
    
    [~,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',sensing_dir',ch_types,Lin); %calculating SSS matrix
    conditionNumber = cond(SNin);

else
    sensing_dir = cos(angle).*phi_hat + sin(angle).*theta_hat;
    % phi_hat(:,2)= angle;
    
    [~,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix',R_hat',theta_hat',sensing_dir',ch_types,Lin); %calculating SSS matrix
    conditionNumber = cond(SNin);
end