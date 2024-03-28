% Convert PQW states to Cartesian states
function [position_vec_, velocity_vec_, spef_enegry_] = orbitalElementsToCart(a_, e_, i_param_, omega_param_, w_param_, true_anom_, mew_, deg_or_rad)
    
    if deg_or_rad == "deg"
        i_param_ = i_param_*pi/180;
        omega_param_ = omega_param_*pi/180;
        w_param_ = w_param_*pi/180;
        true_anom_ = true_anom_*pi/180;        
    end
    % find r_pqw_vec
    r = (a_*(1-e_^2))/(1+e_*cos(true_anom_)); % First get current radius at true anomaly
    r_pqw_vec = [r*cos(true_anom_);r*sin(true_anom_);0]; % Now get Rpqw vector
    
    % Find semi-latus rectum (P)
    p = a_*(1-e_^2);
    
    % Get Vpqw
    v_pqw_vec = [(sqrt(mew_/p)*-sin(true_anom_));(sqrt(mew_/p)*(e_+cos(true_anom_)));0];

    % Get transformation matrix
    Rz_omega = [cos(-omega_param_), sin(-omega_param_),0; sin(-omega_param_), cos(-omega_param_),0; 0,0,1];
    Rx_i = [1,0,0; 0, cos(-w_param_), sin(-w_param_); 0, sin(-w_param_), cos(-w_param_)];
    Rz_w = [cos(-w_param_), sin(-w_param_), 0 ;sin(-w_param_), cos(-w_param_), 0; 0,0,1];

    Tpqw_cart = Rz_omega*Rx_i*Rz_w;
    r_cart_vec = Tpqw_cart*r_pqw_vec;
    v_cart_vec = Tpqw_cart*v_pqw_vec;

    position_vec_ = r_cart_vec;
    velocity_vec_ = v_cart_vec;
    spef_enegry_ = -mew_/(2*a_);
    
end
