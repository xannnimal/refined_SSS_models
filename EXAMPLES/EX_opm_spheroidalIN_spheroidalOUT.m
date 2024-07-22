%% EXAMPLE: OPM sensors, spheroid-in, spheroid-out

%% constant variables regardless
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% opm geometry from Peter at SANDIA
filename="headwithsensors1.mat";
%generate helmet pos and ori with "gen_opm_geometry"
[opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename);
nchan = size(ch_types,1);

%% SSS expansions- multi origin interior
%speficy sensing direction. SQUID=R_hat or EZ, OPM=Theta or phi hat
sensing_dir=R_hat;
other_dir2 = phi_hat;
other_dir1 = theta_hat;
%find major and minor axis of spheroidal ellipse
%Y must be the longest axis of opm_matrix
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_matrix);
%calculate spheroidal in and single-vsh out
%matrix,chan_ori,semi_major,semi_minor,Lin,Lout
[Sin_spm,Sout_spm] = spheroidIN_spheroidOUT(opm_matrix,sensing_dir,origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm,2)
  SNin_spm(:,j) = Sin_spm(:,j)/norm(Sin_spm(:,j));
end

for j = 1:size(Sout_spm,2)
  SNout_spm(:,j) = Sout_spm(:,j)/norm(Sout_spm(:,j));
end

%% generate single dipole simulated data
%current dipole using Samu's implementation of Sarvas
mags = 1:1:nchan;
rs=[0,0,0];
q=[0,1,0]; %y direction
r0=[0.05,0,0]; %5cm along x axis
%phi_0p = dipole_field_sarvas(rs',q',r0',opm_matrix',R_hat',theta_hat',phi_hat',mags)';

dip_mom_out=[1,0,0];
dip_pos_out = [0,0,1.5];
f_start = 100; % start frequency
f_end = 50; % end frequency
f_start_out = 50; % start frequency
f_end_out = 30; % end frequency
timestep = 0.0001;
T = 0.05;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;

for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2));
end

%current dipole in, magnetic dipole out
for i=(1:size(times,2))
    %phi_in_c(:,i) = current_dipole(R',EX',EY',EZ',dip_pos, dip_mom_t(:,i), ch_types)';
    %phi_in(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos',dip_mom_t(:,i),ch_types)'; 
    phi_out(:,i) = magneticDipole(opm_matrix',other_dir1',other_dir2',sensing_dir',dip_pos_out', dip_mom_t_out(:,i),ch_types)'*1e9;
    phi_out_point(:,i) = magneticDipole_pointMags(opm_matrix',sensing_dir',dip_pos_out', dip_mom_t_out(:,i))'*1e9;
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',opm_matrix',other_dir1',other_dir2',sensing_dir',mags)';
    phi_in_point(:,i) = dipole_field_sarvas_pointmags(rs',q_t(:,i),r0',opm_matrix',other_dir1',other_dir2',sensing_dir')';
end 

rng(1,'twister')
noise = randn(size(phi_in,1),size(phi_in,2));
amplitude = 0.15 * phi_in;
%%%% modify this line to do only internal, in+ext, or add noise %%%
phi_0 = (phi_in + phi_out)*100; %+ amplitude .* noise
phi_0point = phi_in_point + phi_out_point;

%% reconstrct internal data
pS=pinv([SNin_spm SNout_spm]);   
XN=pS*phi_0;
data_rec=real(SNin_spm*XN(1:size(SNin_spm,2),:)); 
XN_point=pS*phi_0point;
data_rec_point=real(SNin_spm*XN_point(1:size(SNin_spm,2),:)); 
%check condition numbers
condition_in = cond(SNin_spm);
condition_out= cond(SNout_spm);
condition_both = cond([SNin_spm SNout_spm]);

%check subspace angle
oid_oid=[SNin_spm SNout_spm];
for i=(1:size(times,2))
    check_data_oid_oid(i) = subspace(phi_0(:,i), oid_oid)*180/pi;
    check_data_oid_oid_point(i) = subspace(phi_0point(:,i), oid_oid)*180/pi;
end
check_data_oid_oid_min = min(check_data_oid_oid);
check_data_oid_oid_max = max(check_data_oid_oid);
check_data_oid_oid_av = mean(check_data_oid_oid);

check_data_oid_oid_min_p = min(check_data_oid_oid_point);
check_data_oid_oid_max_p = max(check_data_oid_oid_point);
check_data_oid_oid_av_p = mean(check_data_oid_oid_point);

%% plot data to check
%plot data from single channel
% chan_num=1; %3 corresponds to MEG0111
% data_time=dipole_data.time{1,1};
% data_chan_num=dipole_data.trial{1,1}(chan_num,:); 
% 
% figure(2);
% hold on;
% plot(data_time(:,1:100), data_chan_num(:,1:100))
% plot(data_time(:,1:100),data_rec(chan_num,1:100))
% title('spheroidal and vsh, Channel 1, Sandia Helmet phi, dipole 1cm x')
% xlabel('time')
% ylabel('MEG0121')
% %ylim([-8e-12 8e-12])
% legend({'Raw Data','Reconstructed'},'location','northwest')
% hold off

