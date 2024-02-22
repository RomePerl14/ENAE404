
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

    % Create unit vectors that we will utilize later
    i_hat = [1;0;0];
    j_hat = [0;1;0];
    k_hat = [0;0;1];
    
    % Get the specific angular momentum
    h = cross(position_vector_,velocity_vector_);
    norm_h = norm(h);
    
    % Get the node vector 
    n = cross(k_hat, h);
    norm_n = norm(n);
    n_hat = n/norm_n;
    
    % Get the eccentricity
    e = (cross(velocity_vector_,h)/mew_) - (position_vector_/norm(position_vector_));
    norm_e = norm(e);
    e_x_ = e(1);
    e_y_ = e(2);
    e_z_ = e(3);
    
    % Get the position and velocity vectors magnitudes
    norm_vel = norm(velocity_vector_);
    norm_pos = norm(position_vector_);
    
    % Get the specific energy
    spef_energy_ = (norm_vel^2)/2 - mew_/norm(position_vector_);
    
    % Get a param
    a_ = (-mew_)/(2*spef_energy_);

    % Get i param
    i_param_ = acos((dot(h,k_hat))/norm_h);
    
    % Get omega param
    omega_param_ = acos(dot(n, i_hat)/norm_n);
    
    % Get w param
    w_param_ = acos(dot(n, e)/(norm_n*norm_e));
    
    % Get the true anomaly
    true_anom_ = acos(dot(e,position_vector_)/(norm_e*norm_pos));


end
