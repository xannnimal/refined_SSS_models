function [opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry_avg_full(filename)
%% given a matrix (x,y,z) of OPM sensors, find r-hat, phi-hat, and theta-hat
% Xan McPherson, 2023
%INPUT filename of SANDIA helmet geometry, ".mat" file
%OUTPUT r_hat,theta_hat,phi_hat three nchanx3 matricies normal coil ori

geom = load(filename);
opm_geometry = geom.Headwithsensors1ScaledHeadPoints; %table

opm_matrix=[opm_geometry.x,opm_geometry.y,opm_geometry.z]; %matrix
opm_matrix = order_opm_data(opm_matrix);
opm_avg_len = length(opm_matrix);
opm_matrix_avg = zeros(opm_avg_len,3);

arr = zeros(opm_avg_len,10);
RX = arr(1);
RY= arr(2);
RZ= arr(3);
theta= arr(4);
theta_x= arr(5);
theta_y= arr(6);
theta_z= arr(7);
phi= arr(8);
phi_x= arr(9);
phi_y= arr(10);

for j=1:(opm_avg_len/4)
    if j==1
        avg = mean(opm_matrix(j:4*j,:));
        opm_matrix_avg(j:4*j,:) = repmat(avg,4,1);
    else
        avg = mean(opm_matrix(4*j-3:4*j,:));
        opm_matrix_avg(4*j-3:4*j,:) = repmat(avg,4,1);
    end
end

xsens=opm_matrix_avg(:,1); 
ysens=opm_matrix_avg(:,2); 
zsens=opm_matrix_avg(:,3); 

% find r-hat vector
for i=(1:size(xsens,1))
    RX(i)=xsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
    RY(i)=ysens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
    RZ(i)=zsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
end
R_hat=[RX',RY',RZ']; 

%find theta and phi hat
for i=(1:size(xsens,1))
    theta(i)=acos(zsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2));
    phi(i)= sign(ysens(i))*acos(xsens(i)/sqrt(xsens(i)^2+ ysens(i)^2));
end
phi_z=zeros(1,size(xsens,1));
for i=(1:size(xsens,1))
    theta_x(i)=cos(theta(i))*cos(phi(i));
    theta_y(i)=cos(theta(i))*sin(phi(i));
    theta_z(i)=-sin(theta(i));

    phi_x(i)= -sin(phi(i));
    phi_y(i)=cos(phi(i));
end
theta_hat=[theta_x',theta_y',theta_z'];
phi_hat=[phi_x',phi_y',phi_z'];

ch_types = ones(size(xsens,1),1); % 1 = magnetometers, 0 = gradiometers for S_in code

end