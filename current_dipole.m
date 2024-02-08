function magneticFlux =current_dipole(chanpos,EX,EY,EZ,dip_pos, Q, ch_types)
% Xan 2024
% Xan 2024
% Calculating the magnetic field from a dipole using J Sarvas eq. 25
% current dipole Q at r_o (dip_pos)
% chanpos: nchanx3 of sensor positions
% EX,EY,EZ: nchanx3 each, coil orientation matricies
% dip_pos: 1x3 position vector, from the center of the magnetic dipole to where the magnetic field is measured
% Q: 1x3 current dipole
ch_types=ch_types';
mag_size= 21.0e-3;
weights_mag = [16/81;25/324;25/324;25/324;25/324;10/81;10/81;10/81;10/81];
d = sqrt(3/5)*(mag_size/2);
D_mag = [0 0; d d; -d d; -d -d; d -d; 0 d; 0 -d; d 0; -d 0]';
baseline = 16.69e-3;
weights_grad = (1/(4*baseline))*ones(4,1);
weights_grad(5:8) = -(1/(4*baseline))*ones(4,1);
dx1 = 5.89e-3;
dx2 = 10.8e-3;
dy = 6.71e-3;
D_grad = [dx1 dy; dx2 dy; dx1 -dy; dx2 -dy; -dx1 dy; -dx2 dy; -dx1 ...
   -dy; -dx2 -dy]';

%rhat = dip_pos/norm(dip_pos);
nchan=size(chanpos,1);
for k=(1:nchan)
    if ch_types(k)==1 %magnetometer calculations
        %expressing the location vector r as the product of its magnitude times the unit vector in its direction
        %rhat points from dipole location to sensor
        magneticFlux(k) = mag(chanpos(k,:),EX(k,:),EY(k,:),EZ(k,:),dip_pos, Q,weights_mag, D_mag);
    else %gradiometers, ch_types=0
        magneticFlux(k) = grad(chanpos(k,:),EX(k,:),EY(k,:),EZ(k,:),dip_pos, Q,weights_grad, D_grad);
    end
end


%calculate grads
%In the current dipole calculation, you first take the dot product between (Q x r_0) and r, then multiply del-F by this scalar. And yes, Q, is a single 3D vector just like m.
function mag_flux = grad(chanpos_i,EX_i,EY_i,EZ_i,dip_pos, Q,weights_grad, D_grad)
    mu0 = pi*4E-7;  
    for i=1:8
        r_integral = chanpos_i + D_grad(1,i)'*EX_i + D_grad(2,i)'*EY_i;
        r_vec = r_integral - dip_pos;
        a_vec(i,:)=bsxfun(@minus, r_vec, dip_pos); %r(i,:)-dip_pos;
        a=norm(a_vec); %magnitude
        r=norm(r_vec);
        F = a*(r*a+r^2-dot(dip_pos,r_vec));
        del_F=((1/r)*a^2+(1/a)*dot(a_vec(i,:),r_vec)+2*a+2*r)*r_vec - (a+2*r+(1/a)*dot(a_vec(i,:),r_vec)*dip_pos);
        Bfield(i,:)= (mu0/(4*pi*F^2))*(cross(F*Q,dip_pos)-dot(cross(Q,dip_pos),r_vec)*del_F);
        %Bfield2= Q/r_hat^3;
        B(i,:)= Bfield(i,:);
    end
    mag_flux = dot(B'*weights_grad,EZ_i');
end

%calculate mags
function mag_flux = mag(chanpos_i,EX_i,EY_i,EZ_i,dip_pos, Q,weights_mag, D_mag)
    mu0 = pi*4E-7; 
    for i=1:9
        r_integral = chanpos_i + D_mag(1,i)'*EX_i + D_mag(2,i)'*EY_i;
        r_vec = r_integral - dip_pos;
        a_vec(i,:)=bsxfun(@minus, r_vec, dip_pos); %r(i,:)-dip_pos;
        a=norm(a_vec); %magnitude
        r=norm(r_vec);
        F = a*(r*a+r^2-dot(dip_pos,r_vec));
        del_F=((1/r)*a^2+(1/a)*dot(a_vec(i,:),r_vec)+2*a+2*r)*r_vec - (a+2*r+(1/a)*dot(a_vec(i,:),r_vec)*dip_pos);
        Bfield(i,:)= (mu0/(4*pi*F^2))*(cross(F*Q,dip_pos)-dot(cross(Q,dip_pos),r_vec)*del_F);
        %Bfield2= Q/r_hat^3;
        B(i,:)= Bfield(i,:);
    end
    mag_flux = dot(B'*weights_mag,EZ_i');
end

end