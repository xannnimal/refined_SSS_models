function magneticField = magneticDipole(r,m)
%Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% m: magnetic moment
% r: position vector, from the center of the magnetic dipole to where the magnetic field is measured
mu0 = pi*4E-7; 
mag_r = norm(r); % magnitude of r vector
rhat = r/mag_r; % unit vector of r vector
magneticField = (mu0/(4*pi*mag_r^3))*(3*(dot(m,rhat)*rhat-m));

end
