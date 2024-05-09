
%% Q2
% Find stuff with things
mew_mars = 0.042828*10^6;

% Given
rp = 1000;
e = 0.25;
a = rp/(1-e);
ra = a*(1+e);

% Get the periapsis speed
v_peri = sqrt((2*mew_mars)/rp - (2*mew_mars)/(rp+ra));
v_circ = sqrt(mew_mars/rp);
deltaV = v_peri-v_circ;

% now solve for propellent change
Isp = 250;
g = 9.81*(10^-3);
mi = 1500;
mi_over_mf = exp(deltaV/(Isp*g));
mf = mi/mi_over_mf;

deltaM = mi*(1-exp(-deltaV/(Isp*g)));
mf-mi;

%%% part 2
% given:
struct_ratio = 0.15;

% Find n (mass ratio)
n = mi/mf;

% now, find lambda (payload ratio)
lambda = (n*struct_ratio-1)/(1-n);

% For shits and giggles, get payload mass
mpl = (lambda*mi)/(1+lambda);
me = struct_ratio*(mi-mpl);

