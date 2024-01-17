function B_field = current_dipole(chanpos, dip_pos, Q)
% Iman and Xan 2024
%Calculating the magnetic field from a dipole using J Sarvas eq. 25
% current dipole Q at r_o (dip_pos)
% dip_pos: position vector, from the center of the magnetic dipole to where the magnetic field is measured
% chanpos: channel position matrix nchanx3
mu0 = pi*4E-7; 
rhat = dip_pos/norm(dip_pos);
%make r vector from (x,y,z) pos
for i=(1:size(chanpos,1))
    r(i)=sqrt(chanpos(i,1)^2+chanpos(i,2)^2+chanpos(i,3)^2); %1x144, r x channel
end
%calculate F, del_F, and B
for i=(1:size(chanpos,1))
    r_vec=chanpos(i,:);
    c=bsxfun(@minus, r_vec, dip_pos);
    a_vec(i,:)=bsxfun(@minus, r_vec, dip_pos); %r(i,:)-dip_pos;
    a=norm(a_vec); %magnitude
    r=norm(r_vec);
    F = a*(r*a+r^2-dot(dip_pos,r_vec));
    del_F=((1/r)*a^2+(1/a)*dot(a_vec(i,:),r_vec)+2*a+2*r)*r_vec - (a+2*r+(1/a)*dot(a_vec(i,:),r_vec)*r_vec);
    %B_field(i,:)= (mu0/(4*pi*F^2))*(cross(F*Q,dip_pos)-dot(cross(Q,dip_pos),r_vec*del_F));
end


end