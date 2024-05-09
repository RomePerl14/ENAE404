%-------- HW 8 MATLAB code --------%
% Romeo Perlstein, section 0101    %
% Chat it's over... I'm cooked!!!  %

%% Q1
% given the 3-2-1 euler rotation of a body:
theta1= 30;
theta2 = 40;
theta3 = 10;
Tz = [cosd(theta1), -sind(theta1), 0; sind(theta1), cosd(theta1), 0; 0,0,1];
Ty = [cosd(theta2), 0, sind(theta2); 0, 1, 0; -sind(theta2), 0, cosd(theta2)];
Tx = [1,0,0; 0, cosd(theta3), -sind(theta3); 0, sind(theta3), cosd(theta3)];

R_full = Tz*Ty*Tx


% Get the principle rotation angle
phi = acosd(.5*(R_full(1,1)+R_full(2,2)+R_full(3,3)-1));

% Get the principle axis
e_vec = 1/(2*sind(phi)) * [R_full(2,3)-R_full(3,2);R_full(3,1)-R_full(1,3);R_full(1,2)-R_full(2,1)]

% Now, find the quaternion values:
q1 = e_vec(1)*sind(phi/2)
q2 = e_vec(2)*sind(phi/2)
q3 = e_vec(3)*sind(phi/2)
q4 = cosd(phi/2)
if((q1^2 + q2^2 + q3^2 + q4^2) > .99999 && (q1^2 + q2^2 + q3^2 + q4^2) < 1.00001)
    fprintf("Quaternions check out!\n\n");
end


