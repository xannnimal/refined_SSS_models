function magneticFlux = magneticDipole(chanpos,sensing_dir,grad_pos_ex,dip_pos, dip_mom, ch_types)
% Iman and Xan 2024
% Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% dip_mom: magnetic moment 1x3
% dip_pos: location of dipole
% chanpos: channel position matrix nchanx3
% sensing_dir: the orientation/measuring direction of the sensors
% grad_pos_ex: for gradiometers, need EX positions, for OPM sensors this
%                part of the function is never used
ch_types=ch_types';
mu0 = pi*4E-7; 
%rhat = dip_pos/norm(dip_pos);
%make r vector from (x,y,z) chanpos 
for i=(1:size(chanpos,1))
    %r(i)=sqrt(chanpos(i,1)^2+chanpos(i,2)^2+chanpos(i,3)^2);
    sensing_dir(i,:)=sensing_dir(i,:)/norm(sensing_dir(i,:));
end
%calculate B for each channel position
for i=(1:size(chanpos,1))
    if ch_types(i)==1 %magnetometer calculations
        %expressing the location vector r as the product of its magnitude times the unit vector in its direction
        %rhat points from dipole location to sensor
        r_vec=[chanpos(i,1)-dip_pos(1),chanpos(i,2)-dip_pos(1),chanpos(i,3)-dip_pos(3)];
        r= norm(r_vec);
        rhat=[r_vec(1)/r, r_vec(2)/r, r_vec(3)/r];
        BField = mu0/(4*pi)*(1/(r^3))*(3*(dot(dip_mom,rhat)*rhat-dip_mom));
        magneticFlux(i,:) =dot(BField,sensing_dir(i,:)); %no integral for point mags
    else %gradiometers, ch_types=0
        %calculate magntiudes of r
        del_x = 0.0084*grad_pos_ex(i,:);
        plus_r_vec = [chanpos(i,1)+del_x(1)-dip_pos(1), chanpos(i,2)+del_x(2)-dip_pos(2), chanpos(i,3)+del_x(3)-dip_pos(3)];
        plus_r = norm(plus_r_vec);
        plus_rhat = [plus_r_vec(1)/plus_r, plus_r_vec(2)/plus_r, plus_r_vec(3)/plus_r];
        min_r_vec = [chanpos(i,1)-del_x(1)-dip_pos(1), chanpos(i,2)-del_x(2)-dip_pos(2), chanpos(i,3)-del_x(3)-dip_pos(3)];
        min_r = norm(min_r_vec);
        min_rhat = [min_r_vec(1)/min_r, min_r_vec(2)/min_r, min_r_vec(3)/min_r];

        %B field 8.4 mm in the sensor's positive x-direction (EX)
        BField_plus = mu0/(4*pi)*(1/(plus_r^3))*(3*(dot(dip_mom,plus_rhat)*plus_rhat-dip_mom));
        %B field 8.4 mm in the sensor's negative x-direction (EX)
        BField_min = mu0/(4*pi)*(1/(min_r^3))*(3*(dot(dip_mom,min_rhat)*min_rhat-dip_mom));
        %flux = positive - negative, divide by the gradiometer baseline (16.8 mm)
        flux = dot(BField_plus-BField_min,sensing_dir(i,:));
        magneticFlux(i,:)=flux/0.0168;
    end
end

end