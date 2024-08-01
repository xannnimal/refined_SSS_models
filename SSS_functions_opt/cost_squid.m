function conditionNumber = cost_squid(angle,R,EX,EY,ch_types,Lin,r,phi)
%%  input: 102 angles, SSS parameters, spherical coordinate parameters
%%  output: condition number
% fixed sensor positions

%find x, y, z coords from spherical coords
x_sph = r.*sin(angle).*cos(phi);
y_sph = r.*sin(angle).*sin(phi);
z_ph = r.*cos(angle);

EZ = [x_sph,y_sph,z_ph]; %new sensing direction based on angle

[~,SNin]= Sin_vsh_vv([0,0,0]',R',EX',EY',EZ',ch_types,Lin);  %calculate SSS matrix
conditionNumber = cond(SNin); %calculate condition number

end