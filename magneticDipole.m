function magneticField = magneticDipole(chanpos,sensing_dir,grad_pos_ex,dip_pos, dip_mom, ch_types)
% Iman and Xan 2024
%Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% dip_mom: magnetic moment 1x3
% dip_pos: position vector, from the center of the magnetic dipole to where the magnetic field is measured
% chanpos: channel position matrix nchanx3
% sensing_dir: the orientation/measuring direction of the sensors
%grad_pos_ex: for gradiometers, need EX positions
ch_types=ch_types';
mu0 = pi*4E-7; 
rhat = dip_pos/norm(dip_pos);
%make r vector from (x,y,z) pos and coil_hat for rhat of sensing direction
for i=(1:size(chanpos,1))
    r(i)=sqrt(chanpos(i,1)^2+chanpos(i,2)^2+chanpos(i,3)^2);
    %coil_hat(i)=sqrt(coilori(i,1)^2+coilori(i,2)^2+coilori(i,3)^2);
end
%calculate B for each channel position
for i=(1:size(chanpos,1))
    if ch_types(i)==1 %magnetometer calculations
        BField(i,:) = mu0/(4*pi)*(1/(r(i)^3))*(3*(dot(dip_mom,rhat)*rhat-dip_mom));
    else %gradiometers, ch_types=0
        plus_x = 0.0084*grad_pos_ex(i,:);
        min_x = -0.0084*grad_pos_ex(i,:);
        plus_r = sqrt((chanpos(i,1)+plus_x(1,1))^2+(chanpos(i,2)+plus_x(1,2))^2+(chanpos(i,3)+plus_x(1,3))^2);
        min_r = sqrt((chanpos(i,1)+min_x(1,1))^2+(chanpos(i,2)+min_x(1,2))^2+(chanpos(i,3)+min_x(1,3))^2);
        %B field 8.4 mm in the sensor's positive x-direction (EX)
        BField_plus(i,:) = mu0/(4*pi)*(1/(plus_r^3))*(3*(dot(dip_mom,rhat)*rhat-dip_mom));
        %B field 8.4 mm in the sensor's negative x-direction (EX)
        BField_min(i,:) = mu0/(4*pi)*(1/(min_r^3))*(3*(dot(dip_mom,rhat)*rhat-dip_mom));
        %B = positive - negative, divide by the gradiometer baseline (16.8 mm)
        BField(i,:)=(BField_plus(i,:)-BField_min(i,:))/0.0168;
    end
end

%project the Bfield into sensor direction, coilpos rhat
for i=(1:size(chanpos,1))
    magneticField(i,:) =dot(sensing_dir(i,:),BField(i,:));
end
end