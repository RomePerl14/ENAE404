%-------- HW 5 MATLAB code --------%
% Romeo Perlstein, section 0101 %

% IT'S A NEW WEEK - HOPEFULLY I CAN KEEP UP THE WORK AND NOT FALL BEHIND!

%% Q1
% Find deltaV for Mercury to Jupiter transfer using conic sections

% Givens:
mew_mercury = 22031.868551; % km^3/s^2 - FROM JPL
mew_saturn = 126712764.1; % km^3/s^2 - FROM JPL
radius_planet_mercury = 2440.5; % km - from NASA fact sheet
r_craft_mercury = 400+radius_planet_mercury; % km
r_mercury = 57.909*10^6; % semi-major axis - from NASA fact sheet

radius_planet_saturn = 60268; % km - from NASA fact sheet
r_craft_saturn = 10000+radius_planet_saturn; % km
r_saturn = 1432.041*10^6; % semi-major axis - from NASA fact sheet

% since we are assuming circular orbits, we need to find the orbit
% velocity of both planets!
mew_sun = 132712*10^6; % from NASA fact sheet
v_mercury = sqrt(mew_sun/r_mercury);
v_saturn = sqrt(mew_sun/r_saturn);

% find velocities of orbits of spacecraft
v_initial_craft_mercury = sqrt(mew_mercury/r_craft_mercury);
v_final_craft_saturn = sqrt(mew_saturn/r_craft_saturn);

% Find the velocity to transfer from mercury to saturn
v_transfer_peri = sqrt(2*((mew_sun/r_mercury) - (mew_sun/(r_mercury+r_saturn))));
v_transfer_apo = sqrt(2*((mew_sun/r_saturn) - (mew_sun/(r_mercury+r_saturn))));

% Now, lets get the escape velocity from mercury
v_escape_mercury = v_transfer_peri - v_mercury;
v_hyperbola_peri_mercury = sqrt(2*((mew_mercury/r_craft_mercury) + ((v_escape_mercury^2)/2)));

% Now we can find the delta V to get to Saturn
deltaV1 = v_hyperbola_peri_mercury - v_initial_craft_mercury;

% Now do saturn
v_escape_saturn = v_saturn - v_transfer_apo;

% MAKING ASSUMPTION THAT HYPERBOLA PERIAPSIS IS SAME AS PARKING ORBIT
% PERIAPSIS, A. BECAUSE THE PROBLEM DOESN'T SAY WE CAN'T AND B. BECAUSE I
% WOULD NOT BE ABLE TO SUBMIT THE HW ON TIME
v_hyperbola_peri_saturn = sqrt(2*((mew_saturn/r_craft_saturn) + ((v_escape_saturn^2)/2)));
deltaV2 = v_hyperbola_peri_saturn - v_final_craft_saturn;

% now get the supplementary info
mass_mercury = .3301*10^24; % kg
mass_saturn = 568.32*10^24; % kg
mass_sun = 1998500*10^24; % kg
SOI_mercury = r_mercury*(mass_mercury/mass_sun)^(2/5); % km - Matches with Wikipedia!
SOI_saturn = r_saturn*(mass_saturn/mass_sun)^(2/5); % km - Matches with Wikipedia!
deltaV_total = deltaV1+deltaV2;

% Find TOF, assuming ellipse:
a = (r_mercury+r_saturn)/2;
e = 1-(r_mercury/a);
E = pi;
t = sqrt((a^3)/mew_sun)*(pi-0.9223*sin(pi)); % seconds!!!

%% Q2
% do a bunch of stuff I don't have time to finish :/
% given:
e2 = 1.2; 
rp2 = 5380;
mew_mars = 0.042828 *10^6; % km^3/s^2 - from NASA fact sheet
% assuming circular orbit
r_mars = 227.956 * 10^6;
v_mars = sqrt(mew_sun/r_mars);

% find escape velocity and velocity at periapsis
v_escape_mars2 = sqrt(((e2-1)*mew_mars)/rp2); % from notes? I guess it works
v_hyperbola_peri_mars2 = sqrt(2*((mew_mars/rp2) + ((v_escape_mars2^2)/2))); % unused

% Get the hyperbola angle
Beta = acosd(1/e2); % get Beta angle

%%% a
% find v_initial, v_final, and then find deltaV
v_mars_vecA = [0;v_mars]; % Only acts in Y direction
v_escape_mars2_vec = [v_escape_mars2*sind(Beta);v_escape_mars2*cosd(Beta)]; % Get the vector components
v_initalA = v_mars_vecA + v_escape_mars2_vec; % Get the vector sum

v_escape_mars2_vec = [v_escape_mars2*sind(Beta);-v_escape_mars2*cosd(Beta)]; % Get the new escape vector
v_finalA = v_mars_vecA + v_escape_mars2_vec; % Get the vector sum

turn_angle = (180-(2*Beta))*(pi/180); % get the turn angle of the parabola
deltaV_A = sqrt(2*(v_escape_mars2^2)*cos(1-turn_angle)); % get the deltaV using the deltaV eq. from notes

% Now get the norms
v_initalA_mag = norm(v_initalA);
v_finalA_mag = norm(v_finalA);

%%% b
% find v_initial, v_final, and then find deltaV
mars_vec_angle = 30; % Given angle away from velocity vector

% Get velocity vectors
v_mars_vecB = [v_mars*sind(mars_vec_angle);v_mars*cosd(mars_vec_angle)];
v_escape_mars2_vec = [v_escape_mars2*sind(Beta);v_escape_mars2*cosd(Beta)];
v_initalB = v_mars_vecB + v_escape_mars2_vec; % Get vector sum

% get final stuff
v_escape_mars2_vec = [v_escape_mars2*sind(Beta);-v_escape_mars2*cosd(Beta)];
v_finalB = v_mars_vecB + v_escape_mars2_vec; % Get vector sum

% Get deltaV (turn angle does not change!)
deltaV_B = sqrt(2*(v_escape_mars2^2)*cos(1-turn_angle));

% Get magnitudes
v_initalB_mag = norm(v_initalB);
v_finalB_mag = norm(v_finalB);

%%% c
% do pretty much exactly the same thing as above, again
v_mars_vecC = [0; -v_mars];
v_escape_mars2_vec = [v_escape_mars2*sind(Beta);v_escape_mars2*cosd(Beta)];
v_initalC = v_mars_vecC + v_escape_mars2_vec;

v_escape_mars2_vec = [v_escape_mars2*sind(Beta);-v_escape_mars2*cosd(Beta)];
v_finalC = v_mars_vecC + v_escape_mars2_vec;

deltaV_C = sqrt(2*(v_escape_mars2^2)*cos(1-turn_angle));

v_initalC_mag = norm(v_initalC);
v_finalC_mag = norm(v_finalC);


%% Q3


