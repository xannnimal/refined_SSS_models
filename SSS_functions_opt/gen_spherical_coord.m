function [cartesian_matrix,r,theta,phi] = gen_spherical_coord(data)
%% given a matrix (x,y,z) of SQUID magnetometer sensors, find r, phi, theta

if size(data,2) > size(data,1)
    data = data';
end

xdir = data(:,1);
ydir = data(:,2);
zdir = data(:,3);
cartesian_matrix=[xdir,ydir,zdir];

r = sqrt(sum(data.^2,2)); 
theta= acos(zdir./r);
phi= sign(ydir).*acos(xdir./sqrt(xdir.^2+ ydir.^2));

end