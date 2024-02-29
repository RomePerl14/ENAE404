% Convert Cartesian state to PQW state values
function [i_param_, omega_param_, w_param_, true_anom_, e_x_, e_y_, e_z_, a_, spef_energy_] = cartToOrbitalElements(position_vector_, velocity_vector_, mew_)
    format long
    if(isempty(position_vector_))
        fprintf("How tf you do that? You need to input a vector of length 3")
    elif(length(position_vector_) == 1)
        fprintf("Inputted position is not a vector! what the frick man?")
    end
    if(isempty(velocity_vector_))
        fprintf("How tf you do that? You need to input a vector of length 3")
    elif(length(velocity_vector_) == 1)
        fprintf("Inputted position is not a vector! what the frick man?")
    end

    %%%------ Initial Parameter Collection ------%%%
    % Create unit vectors that we will utilize later
    i_hat = [1;0;0];
    j_hat = [0;1;0];
    k_hat = [0;0;1];
    
    % Get the specific angular momentum
    h = cross(position_vector_,velocity_vector_); % specific angular momentum vector
    norm_h = norm(h); % magnitude of specific angular momentum vector
    hx = h(1);
    hy = h(2);
    hz = h(3);
    
    % Get the node vector 
    n = cross(k_hat, h); % node vector
    norm_n = norm(n); % magnitude of the node vector
    n_hat = n/norm_n; % unit vector of the node vector
    nx = n(1); % n vector x
    ny = n(2); % n vector y
    nx = n(3); % n vector z
    
    % Get the eccentricity
    e = (cross(velocity_vector_,h)/mew_) - (position_vector_/norm(position_vector_)); % eccentricity vector
    norm_e = norm(e); % magnitude of eccentricity vector
    e_x_ = e(1); % x value for checking
    e_y_ = e(2); % y value for checking
    e_z_ = e(3); % z value for checking
    
    % Get the position and velocity vectors magnitudes
    norm_vel = norm(velocity_vector_); % velocity magnitude
    norm_pos = norm(position_vector_); % position magnitude

    % get the position and velocity x y z values
    rx = position_vector_(1);
    ry = position_vector_(2);
    rz = position_vector_(3);
    vx = velocity_vector_(1);
    vy = velocity_vector_(2);
    vz = velocity_vector_(3);


    % Get the specific energy
    spef_energy_ = (norm_vel^2)/2 - mew_/norm(position_vector_); % specific energy
    
    % Get apoapsis param
    a_ = (-mew_)/(2*spef_energy_); % apoapsis

    %%%------ Different Cases ------%%%
    % Is this necessary? not really, but does it help? no. But it makes it
    % faster right? no, it probably makes it slower in all honesty. But it
    % will at least be correct, right? yes! (probably not actually)
    if(norm_n == 0 && norm_e == 0) % check if orbit is circular, equitorial
        i_param_ = 0;
        w_param_= 0;
        omega_param_ = 0;
        true_anom_ = acos(rx/norm_r); % get the true anomaly
        if(ry < 0) % check condition
            true_anom_ = (2*pi)-true_anom_; % rewrite the true anomaly to account for shift
        end
    elseif(norm_n == 0 && (norm_e < 1 && norm_e > 0)) % check if orbit is elliptical, equitorial
        i_param_ = 0;
        w_param_ = acos(e_x_/norm_e);
        if(e_y_ < 0) % check condition
            w_param_ = (2*pi)-w_param_; % rewrite the w param to account for shift
        end
        omega_param_ = 0;
        true_anom_ = acos((dot(e, position_vector_))/(norm_e*norm_r));
        if(dot(position_vector_, velocity_vector_) < 0)
            true_anom_ = (2*pi)-true_anom_; % rewrite the true anomaly to account for shift
        end
    elseif(norm_n > 0 && norm_e == 0) % check if orbit is circular inclined
        i_param_ = acos(hz/norm_h);
        w_param_ = 0;
        omega_param_ = acos(nx/norm_n);
        if(ny < 0) % check condition
            omega_param_ = (2*pi)-omega_param_; % rewrite the omega param to account for shift
        end
        true_anom_ = acos((dot(n, position_vector_))/(norm_n*norm_pos));
        if(rz < 0) % check condition
            true_anom_ = (2*pi)-true_anom_; % rewrite the true anomaly to account for shift
        end
    elseif(norm_n > 0 && norm_e > 0) % check if orbit is elliptical inclined
        i_param_ = acos(hz/norm_h);
        w_param_ = acos((dot(n,e))/(norm_n*norm_e));
        if(e_z_ < 0)
            w_param_ = (2*pi)-w_param_; % rewrite w param to account for shift
        end
        omega_param_ = acos(nx/norm_n);
        if(ny < 0)
            omega_param_ = (2*pi)-omega_param_; % rewrite omega param to account for shift
        end
        true_anom_ = acos((dot(e,position_vector_))/(norm_e*norm_pos));
        if(dot(position_vector_,velocity_vector_) < 0)
            true_anom_ = (2*pi)-true_anom_; % rewrite true anomaly to account for shift
        end
    end
end
