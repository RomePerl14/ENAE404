format long
r = [3634.1 ; 5926 ; 1206.6];
v = [-6.9049 ; 4.3136 ; 2.6163];

h = cross(r,v)
norm_h = norm(h)

i_hat = [1;0;0];
j_hat = [0;1;0];
k_hat = [0;0;1];


n = cross(k_hat, h)
norm_n = norm(n)
n_hat = n/norm_n

mew_Earth = 398600.44;

e = (cross(v,h)/mew_Earth) - (r/norm(r))
norm_e = norm(e)

norm_v = norm(v)
norm_r = norm(r)

spef_Energy = (norm_v^2)/2 - mew_Earth/norm(r)

a = (-mew_Earth)/(2*spef_Energy)

i_param = acos((dot(h,k_hat))/norm_h)

omega_param = acos(dot(n, i_hat)/norm_n)

w_param = acos(dot(n, e)/(norm_n*norm_e))

true_anom = acos(dot(e,r)/(norm_e*norm_r))
