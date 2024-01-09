function magneticField = magneticDipole(chanpos,dip_pos,dip_mom)
%Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% dip_mom: magnetic moment 1x3
% dip_pos: position vector, from the center of the magnetic dipole to where the magnetic field is measured
% chanpos: channel position matrix nchanx3
mu0 = pi*4E-7; 
rhat = dip_pos/norm(dip_pos);
%make r vector from (x,y,z) pos
for i=(1:size(chanpos,1))
    r(i)=sqrt(chanpos(i,1)^2+chanpos(i,2)^2+chanpos(i,3)^2);
end
%calculate B for each channel position
for i=(1:size(chanpos,1))
    magneticField_i = mu0/(4*pi)*(1/r(i))*(3*(dot(dip_mom,rhat)*rhat-dip_mom));
    magneticField(i,:)=magneticField_i;
end

% mag_r = norm(r); % magnitude of r vector
% rhat = r/mag_r; % unit vector of r vector
% magneticField = (mu0/(4*pi*mag_r^3))*(3*(dot(m,rhat)*rhat-m));

end
