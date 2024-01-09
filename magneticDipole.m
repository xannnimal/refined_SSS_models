function magneticField = magneticDipole(chan_pos,dip_pos dip_mom)
%Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% m: magnetic moment
% r: position vector, from the center of the magnetic dipole to where the magnetic field is measured
mu0 = pi*4E-7; 
mag_r = norm(dip_pos); % magnitude of r vector
rhat = dip_pos; % unit vector of r vector
chan_pos_r = zeros(height(chan_pos),1);
magneticField = zeros(height(chan_pos),1);
for i = 1:height(chan_pos)
    chan_pos_r(i,1) = sqrt(chan_pos(i,1)^2 + chan_pos(i,2)^2 + chan_pos(i,3)^2);
    magneticField(i,1) = (mu0/(4*pi*chan_pos_r(i,1)^3))*(3*(dot(dip_mom,rhat)*rhat-dip_mom));
    % disp(magneticField(i,1))

magneticField
end
% 144x3 matrix convert x y z into r vector 144x1
