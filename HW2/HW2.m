%-------- HW 2 MATLAB code --------%
% Romeo Perlstein, section 0101 %
%% Q1
% See supporting writeup for Q1

%% Q2
% Use matlab func to find orbital elements

%-- Given values --%
r = [3634.1 ; 5926 ; 1206.6]; % position vec
v = [-6.9049 ; 4.3136 ; 2.6163]; % velocity vec
mew_Earth = 398600.44; % gravitational parameter

% Get the orbital elements from the original data, for shits and giggles
[i_param, omega_param, w_param, true_anom, ex, ey, ez, a, spef_energy] = cartToOrbitalElements(r, v, mew_Earth, "rad")

%%%----- Q2-1 -----%%%
% Plot the changing orbital elements as subplots

%--- ODE func values ---%
tall_er_ant = (10^-13); % Tolerance
step_size = 1000; % step size 
max_time = 70000000; % max time (0->max_time)
t = [0:step_size:max_time]; % timestep

% ODE options
ODE_options = odeset("RelTol", tall_er_ant, "AbsTol", tall_er_ant);

% Didymos Orbit information
mew_sun = 1.32712 * (10^11);
didymos_initial_x = -2.39573*10^8;
didymos_initial_y = -2.35661*10^8;
didymos_initial_z = 9.54384*10^6;
didymos_initial_vx = 1.24732*10^1;
didymos_initial_vy = -9.74427*10^0;
didymos_initial_vz = -8.78661*10^-1;
didymos_initial_state = [didymos_initial_x; didymos_initial_y; didymos_initial_z; didymos_initial_vx; didymos_initial_vy; didymos_initial_vz];

% Get the positions and velocities of the entire orbit in cartesian coords
[T,Y] = ode45(@myodefun, t, didymos_initial_state, ODE_options, mew_sun);

for i=1:length(t)
    [i_didymos(i), omega_didymos(i), w_didymos(i), true_anom_didymos(i), ex_didymos(i), ey_didymos(i), ez_didymos(i), a_didymos(i), spef_energy_didymos(i)] = cartToOrbitalElements([Y(i,1);Y(i,2);Y(i,3)], [Y(i,4);Y(i,5);Y(i,6)], mew_sun, "rad");
    e(i) = sqrt(ex_didymos(i)^2 + ey_didymos(i)^2 + ez_didymos(i)^2);
end

%%%--- Plot the the orbital element states ---%%%
close all
tiledlayout(2,4) % make the tiled layout


nexttile % first subplot - i param
plot(t, (i_didymos*(180/pi)))
axis([0, max_time, 0, 6])
ylabel("Inclination: i param (degrees)")
xlabel("Time (seconds)")
title("Inclination of Didymos vs Time")

nexttile % second subplot - w param
plot(t, (w_didymos*(180/pi)))
% axis([0, max_time, 0, 360])
ylabel("\omega (degrees)")
xlabel("Time (seconds)")
title("\omega of Didymos vs Time")

nexttile % third subplot - omega param
plot(t, (omega_didymos*(180/pi)))
% axis([0, max_time, 0, 360])
ylabel("\Omega (degrees)")
xlabel("Time (seconds)")
title("\Omega of Didymos vs Time")

nexttile % fourth subplot - true anomaly param
plot(t, (true_anom_didymos*(180/pi)))
% axis([0, max_time, 0, 360])
ylabel("True Anomaly (degrees)")
xlabel("Time (seconds)")
title("True Anomaly of Didymos vs Time")

nexttile % fifth subplot - a param
plot(t, (a_didymos*(180/pi)))
% axis([0, max_time, 0, 360])
ylabel("Semi-major axis (km)")
xlabel("Time (seconds)")
title("Semi-major axis of Didymos vs Time")

%%%------ Q2-2 ------%%%
% The plots make sense because the orbital elements i, \omega, and \Omega
% are CONSTANT state values in the PQW frame. The element that changes in
% teh PQW frame is the true anomaly, while the i, \omega and \Omega values
% simply represent the angles of the orbit with respect to the interial
% frame. True Anomaly represents the current angle of the spacecraft!. And
% of course, the semi-major axis and eccentricity should not change in a
% table orbit

%% Q3



% From ENAE301
function ydot = myodefun(t, y, mew)
    r_mag = norm(y(1:3));
    ydot(1,1) = y(4);
    ydot(2,1) = y(5);
    ydot(3,1) = y(6);
    ydot(4,1) = (-mew/r_mag^3)*y(1);
    ydot(5,1) = (-mew/r_mag^3)*y(2);
    ydot(6,1) = (-mew/r_mag^3)*y(3);
end
