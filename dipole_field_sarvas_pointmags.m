% Flux through pick-up loops (magnetometers and gradiometers) based on the
% 1987 formula by J. Sarvas
% Modified for point magnetometer calculations
% vuovec = dipole_field_sarvas(rs,q,r0,R,NX,NY,NZ,mags) where
% rs: center of the conducting sphere
% q: dipole moment vector
% r0: location of the dipole in the coordinate system of the conductor
% R: sensor locations
% NX,NY,NZ: direction vectors of the sensors (z-direction is the normal)
% mags: indices of the magnetometer channels (others are assumed
% gradiometers)

function vuovec = dipole_field_sarvas_pointmags(rs,q,r0,R,NX,NY,NZ)

mag_size = 21.0e-3; % Side length of a square magnetometer loop (in meters)
weights_mag = [16/81;25/324;25/324;25/324;25/324;10/81;10/81;10/81;10/81]; % Weights for the discretized surface integral
d = sqrt(3/5)*(mag_size/2);
D_mag = [0 0; d d; -d d; -d -d; d -d; 0 d; 0 -d; d 0; -d 0]';
dx1 = 5.89e-3;
dx2 = 10.8e-3;
dy = 6.71e-3;

nchan = size(R,2); % Number of channels
for i = 1:nchan
   vuovec(i) = magn_vaste(q,r0,R(:,i)-rs,NX(:,i),NY(:,i),NZ(:,i),weights_mag,D_mag);
end


function vuo = magn_vaste(q,r0,rm,nx,ny,nz,weights_mag,D_mag)
rint = rm; %+ D_mag(1,i)*nx + D_mag(2,i)*ny;
a = rint - r0;
na = norm(a);
nr = norm(rint);
F = na*(nr*na + nr^2 - dot(r0,rint));
NABLA_F = ((1/nr)*na^2 + (1/na)*dot(a,rint) + 2*na + 2*nr)*rint - ...
(na + 2*nr + (1/na)*dot(a,rint))*r0;
qr = cross(q,r0);
B= (F*qr - dot(qr,rint)*NABLA_F)/F^2;
vuo = dot(B,nz); %no B*weights_mag





