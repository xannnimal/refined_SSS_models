function Flux = current_dipole_pointmags(chanpos, coilori, dip_pos, Q)
% Xan 2024
% Calculating the magnetic field from a dipole using J Sarvas eq. 25
% current dipole Q at r_o (dip_pos)
% chanpos: nchanx3 of sensor positions
% coilori: nchanx3 sensing direction coil matrix
% dip_pos: 1x3 position vector, from the center of the magnetic dipole to where the magnetic field is measured
% Q: 1x3 current dipole
mu0 = pi*4E-7; 

%calculate F, del_F, and B
%In the current dipole calculation, you first take the dot product between (Q x r_0) and r, then multiply del-F by this scalar. And yes, Q, is a single 3D vector just like m.
for i=(1:size(chanpos,1))
    r_vec=chanpos(i,:);
    a_vec(i,:)=bsxfun(@minus, r_vec, dip_pos); %r(i,:)-dip_pos;
    a=norm(a_vec); %magnitude
    r=norm(r_vec);
    F = a*(r*a+r^2-dot(dip_pos,r_vec));
    del_F=((1/r)*a^2+(1/a)*dot(a_vec(i,:),r_vec)+2*a+2*r)*r_vec - (a+2*r+(1/a)*dot(a_vec(i,:),r_vec)*dip_pos);
    B(i,:)= (mu0/(4*pi*F^2))*(cross(F*Q,dip_pos)-dot(cross(Q,dip_pos),r_vec)*del_F);
end

%project the Bfield into sensor direction, coilpos rhat
for i=(1:size(chanpos,1))
    Flux(i,:) =dot(coilori(i,:),B(i,:));
end

end
