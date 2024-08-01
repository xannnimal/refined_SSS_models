coordsys = 'device'; 
rawfile = "sample_audvis_raw.fif";

%% generate squid geometry 
[R,EX,EY,EZ, ch_types] = gen_squid_geometry(rawfile,coordsys); 
R = R';
EX = EX';
EY = EY';
EZ = EZ'; % magnetometer sensing direction, normal to helmet
Lin = 8;

[~,r,theta,phi] = gen_spherical_coord(EZ); %find r, theta, phi for initial sensing direction
sensor_len = length(R);

objFun = @(angles) cost_squid(angles,R,EX,EY,ch_types,Lin,r,phi); %objective function
lb = -pi*ones(sensor_len,1); %lower bound
ub = pi*ones(sensor_len,1); %upper bound
x0 = theta; %starting point, theta of initial sensing direction

%run for 1000 iterations, plot condition number
options = optimoptions('simulannealbnd','Display','iter','MaxIterations',10,'PlotFcn',{@saplotbestx,@saplotbestf,@saplotx,@saplotf});
[angles,condition_num_opt] = simulannealbnd(objFun,x0,lb,ub,options);

x_sph = r.*sin(angles).*cos(phi);
y_sph = r.*sin(angles).*sin(phi);
z_sph = r.*cos(angles);
EZ_opt = [x_sph,y_sph,z_sph]';