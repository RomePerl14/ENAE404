%-------- HW 2 MATLAB code --------%
% Romeo Perlstein, section 0101 %

%% Q2
%%% Use matlab func to find orbital elements:

% Given:

r = [3634.1 ; 5926 ; 1206.6];
v = [-6.9049 ; 4.3136 ; 2.6163];
mew_Earth = 398600.44;

[i_param, omega_param, w_param, true_anom, ex, ey, ez, a, spef_energy] = cartToOrbitalElements(r, v, mew_Earth)

%%% Q2-1
% Plot the changing orbital elements as subplots

% Plotting constants
tall_er_ant = (10^-13);
step_size = 10000;
max_time = 70000000;

% Time step
t = [0:step_size:max_time];

% ODE options
ODE_options = odeset("RelTol", tall_er_ant, "AbsTol", tall_er_ant);

% Didymos Orbit information
didymos_initial_x = -2.39573*10^8;
didymos_initial_y = -2.35661*10^8;
didymos_initial_z = 9.54384*10^6;
didymos_initial_vx = 1.24732*10^1;
didymos_initial_vy = -9.74427*10^0;
didymos_initial_vz = -8.78661*10^-1;
didymos_initial_state = [didymos_initial_x; didymos_initial_y; didymos_initial_z; didymos_initial_vx; didymos_initial_vy; didymos_initial_vz;0;0;0];



for i=1:length(t)
    [i_param(i), omega_param(i), w_param(i), true_anom(i), ex(i), ey(i), ez(i), a(i), spef_energy(i)] = cartToOrbitalElements(r, v, mew_Earth);
    e(i) = sqrt(ex(i)^2 + ey(i)^2 + ez(i)^2);
end
tiledlayout(2, 4)
nexttile
plot(0:length(t)-1, i_param)
title("Plot of i Value Over Time")
xlabel("Time (s)")
ylabel("Radians (rads/sec")

nexttile
plot(0:length(t)-1, omega_param)
title("Plot of Omega Value Over Time")
xlabel("Time (s)")
ylabel("Radians (rads/sec")

nexttile
plot(0:length(t)-1, w_param)
title("Plot of Omega Value Over Time")
xlabel("Time (s)")
ylabel("Radians (rads/sec")

nexttile
plot(0:length(t)-1, w_param)
title("Plot of w Value Over Time")
xlabel("Time (s)")
ylabel("Radians (rads/sec")

nexttile
plot(0:length(t)-1, true_anom)
title("Plot of True Anomaly Over Time")
xlabel("Time (s)")
ylabel("Radians (rads/sec")

nexttile
plot(0:length(t)-1, e)
title("Plot of Magnitude of e Over time")
xlabel("Time (s)")
ylabel("Radians (rads/sec")

nexttile
plot(0:length(t)-1, a)
title("Plot of 'a' Over time")
xlabel("Time (s)")
ylabel("")

nexttile
plot(0:length(t)-1, spef_energy)
title("Plot of a spef_energy over time")
xlabel("Time (s)")
ylabel("")
