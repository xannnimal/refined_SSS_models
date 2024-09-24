%% testing spheroid fit
% test spheroidal functions with a spherical sensor array
% should get same results as with spherical SSS
% DO GET SAME RESULTS 9/24/2024
% ADDPATH gensph_equidist is in "invest multiorigin SSS -sensor array creation"
clear

%% generate spherical sensor array
r=0.11; %cm
n= 3;
[xsens,ysens,zsens,d,N_count] = gensph_equidist(r,n);
R = [xsens,ysens,zsens];
% find r-hat vector
for i=(1:size(xsens,1))
    RX(i)=xsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
    RY(i)=ysens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
    RZ(i)=zsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
end
EZ=[RX',RY',RZ']; 

%find theta and phi hat
for i=(1:size(xsens,1))
    theta(i)=acos(zsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2));
    phi(i)=sign(ysens(i))*acos(xsens(i)/sqrt(xsens(i)^2+ ysens(i)^2));
end
phi_z=zeros(1,size(xsens,1));
for i=(1:size(xsens,1))
    theta_x(i)=cos(theta(i))*cos(phi(i));
    theta_y(i)=cos(theta(i))*sin(phi(i));
    theta_z(i)=-sin(theta(i));

    phi_x(i)= -sin(phi(i));
    phi_y(i)=cos(phi(i));
end
EX=[theta_x',theta_y',theta_z'];
EY=[phi_x',phi_y',phi_z'];

ch_types = ones(size(xsens,1),1); % 1 = magnetometers, 0 = gradiometers for S_in code
for i=(1:306)
    mags(i)=i;
end
% check angles and plot sensor coil orientations
% opm_matrix=[xsens,ysens,zsens];
% figure(1)
% hold on;
% grid on;
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),R_hat(:,1),R_hat(:,2),R_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),theta_hat(:,1),theta_hat(:,2),theta_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),phi_hat(:,1),phi_hat(:,2),phi_hat(:,3))
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'r')
% scatter3(0,0,0,'g*')
% %scatter3(center2(1),center2(2),center2(3),'b*')
% title('OPM Helmet')
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% zlabel('z axis (m)')
% hold off;
% rotate3d
% view(135, 20);

%calculate spheroidal in/out
%find semi major and minor
Lin=8;
Lout=3;
[semi_major,semi_minor,origin]=find_ellipse_axis(R);
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R,EZ,origin,semi_major,semi_minor,Lin,Lout);

for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end

%calculate single in/out
R=R';
EX=EX';
EY=EY';
EZ=EZ';
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);

%% generate dependent dipoles
%current dipole using Samu's implementation of Sarvas
rs=[0,0,0];
r0=[0,0,5];
q=[0,1,0]; %y direction
%r0=[0.05,0,0]; %5cm along x axis
%phi_0 = dipole_field_sarvas(rs',q',r0',R,EX,EY,EZ,mags)';

%add time dependence to dipole moment
dip_mom_out=[1,0,0];
dip_pos_out = [0,0,1.5]; %1.5 meters
f_start = 80; % start frequency
f_end = 30; % end frequency
f_start_out = 50; % start frequency
f_end_out = 30; % end frequency
timestep = 0.0001;
T = 0.1;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;
% 
for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2))*5e8;
end
% 
% %current dipole in, magnetic dipole out
for i=(1:size(times,2))
    %phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
    %phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)'; 
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t_out(:,i),ch_types)'*1e-12;
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)'*1e-12;
end
%phi_0=phi_in+phi_out;
%add gaussian noise at 10 percent of max value of phi_0
noise = randn(size(phi_in,1),size(phi_in,2));
amplitude = 0.15 * phi_in;
%%%% modify this line to do only internal, in+ext, or add noise %%%
phi_0 = phi_in + phi_out; % + amplitude .* noise; %
%%%%

%% reconstruct 
%single in, single out
pS=pinv([SNin SNout]);
XN=pS*phi_0;
data_rec_vsh=real(SNin*XN(1:size(SNin,2),:));
%spheroidal in, spheroidal out
pS_sph_sph=pinv([SNin_spm SNout_spm]);   
XN_sph_sph=pS_sph_sph*phi_0;
data_rec_sph_sph=real(SNin_spm*XN_sph_sph(1:size(SNin_spm,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh=pinv([SNin_spm SNout]);   
XN_sph_vsh=pS_sph_vsh*phi_0;
data_rec_sph_vsh=real(SNin_spm*XN_sph_vsh(1:size(SNin_spm,2),:));

%% check condition numbers
cond_vsh_vsh=cond([SNin SNout]);
cond_SNin=cond(SNin);
cond_sph_sph = cond([SNin_spm SNout_spm]);
cond_sph_vsh = cond([SNin_spm SNout]);

%% subsapce angles
sVSH_sVSH=[SNin SNout];
oid_oid=[SNin_spm SNout_spm];
oid_sVSH=[SNin_spm SNout];

for i=(1:size(times,2))
    check_data_vsh_vsh_d(i) = subspace(phi_0(:,i), sVSH_sVSH)*180/pi;
    check_data_oid_oid_d(i) = subspace(phi_0(:,i), oid_oid)*180/pi;
    check_data_oid_vsh_d(i) = subspace(phi_0(:,i), oid_sVSH)*180/pi;
end
check_data_vsh_vsh_dmin = min(check_data_vsh_vsh_d);
check_data_vsh_vsh_dmax = max(check_data_vsh_vsh_d);
check_data_vsh_vsh_dav = mean(check_data_vsh_vsh_d);

check_data_oid_oid_dmin = min(check_data_oid_oid_d);
check_data_oid_oid_dmax = max(check_data_oid_oid_d);
check_data_oid_oid_dav = mean(check_data_oid_oid_d);

check_data_oid_vsh_dmin = min(check_data_oid_vsh_d);
check_data_oid_vsh_dmax = max(check_data_oid_vsh_d);
check_data_oid_vsh_dav = mean(check_data_oid_vsh_d);
