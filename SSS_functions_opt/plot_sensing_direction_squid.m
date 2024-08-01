coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";

[R,EX,EY,EZ, ch_types] = gen_squid_geometry(rawfile,coordsys); 
[cartesian_matrix,r,theta,phi] = gen_spherical_coord(EZ);
angles_opt = load("data/squid_angles_1000.mat");
angles_opt = angles_opt.angles;

x_sph = r.*sin(angles_opt).*cos(phi);
y_sph = r.*sin(angles_opt).*sin(phi);
z_sph = r.*cos(angles_opt);
opt_EZ = [x_sph,y_sph,z_sph]';
subspaceAngle = zeros(length(EZ),1);
for i = 1:length(EZ)
    subspaceAngle(i) = subspace(EZ(:,i),opt_EZ(:,i));
end
subspaceAngle_deg = rad2deg(subspaceAngle);


% figure(2);
% sensors = linspace(1,102,102);
% bar(sensors,subspaceAngle_deg,'EdgeColor','none')


fontSize = 17;
figure(3)
hold on;
histogram(subspaceAngle_deg,'NumBins',10,'FaceColor','#77AC30','FaceAlpha', 0.72);
ax = gca;
ax.FontSize = 15; 
xlabel(['Deviation from initial sensor orientation (' char(176) ')'], 'FontSize', fontSize)
ylabel('Number of Sensors', 'FontSize', fontSize)
hold off;
% 
% fontSize = 17;
% figure(1)
% hold on;
% grid on;
% quiver3(R(1,:),R(2,:),R(3,:),EZ(1,:),EZ(2,:),EZ(3,:),'DisplayName','Initial Sensing Vector')
% quiver3(R(1,:),R(2,:),R(3,:),x_sph',y_sph',z_sph','Color','#77AC30','DisplayName','Optimized Sensing Vector')
% scatter3(R(1,:),R(2,:),R(3,:),'r','DisplayName','Sensor Location')
% ax = gca;
% ax.FontSize = 15; 
% % title('Neuromag 306-Channel Sensor Helmet', 'FontSize', 20)
% xlabel('x-axis (m)', 'FontSize', fontSize)
% ylabel('y-axis (m)', 'FontSize', fontSize)
% zlabel('z-axis (m)', 'FontSize', fontSize)
% legend()
% hold off;
% rotate3d
% view(135, 20);
