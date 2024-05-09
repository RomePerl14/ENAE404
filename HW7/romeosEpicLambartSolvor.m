%--------- Lambert Solver ---------%
%          Romeo Perlstein         %
% Implementation of Prof. Barbee's %
%   Lambert Solver algorithm from  %
%   ENAE404 - Astrodynamics (UMD)  %




function [vel_initial_vec, vel_final_vec] = romeosEpicLambartSolvor(pos_vec_initial_, pos_vec_final_, TOF_, direction_method_, mew_)
    %% SETUP
    % First, get the norm of the inputted position vectors
    pos_vec_initial = pos_vec_initial_;
    pos_vec_final = pos_vec_final_;
    r_initial = norm(pos_vec_initial);
    r_final = norm(pos_vec_final);
    TOF = TOF_;
    mew = mew_;

    % Now, get the cos(deltaV)
    cos_delta_true_anom = (dot(pos_vec_initial,pos_vec_final))/(r_initial*r_final);
    
    % Get the initial A value
    try
        if(direction_method_ == "short")
            A = 1*(sqrt(r_initial*r_final*(1+cos_delta_true_anom)));
        elseif(direction_method_ == "long")
            A = 1*(sqrt(r_initial*r_final*(1+cos_delta_true_anom)));
        else
            fprintf("\nYou need to give a direction method! Either\n'short'\nor\n'long'\nALL LOWERCASE\n\n");
            return
        end
    catch
        fprintf("\nYou need to give a direction method! Either\n'short'\nor\n'long'\nALL LOWERCASE\n\n");
        return
    end

    % Check a few things:
    if(A == 0 || isnan(A))
        fprintf("\nThere is an error! A = 0 or NaN, which should not be!\n\n")
        return
    end
    delta_true_anom = acos(cos_delta_true_anom);
    if(delta_true_anom == 0 || isnan(A))
        fprintf("\nThere is an error! deltaV = 0 or NaN, which should not be!\n\n")
        return
    end

    % If the previous is not the case, then move on
    % Starting values for the loop
    psi = 0;
    C2 = 1/2;
    C3 = 1/6;
    psi_up = 4*(pi^2);
    psi_low = -4*(pi^2);
    deltaT = 0;

    %% Start the while loop to iteratively solve for everything
    tall_er_ant = 0.000001;
    while(abs(TOF-deltaT) >= tall_er_ant)
        % Get y
        y = r_initial + r_final + ( A*((psi*C3)-1) )/sqrt(C2);
        
        % Check our values
        if(A > 0 && y < 0)
            psi_low = psi_low + (pi/4);
        end
        
        % get x and deltaT
        x = sqrt(y/C2);
        deltaT = ( (x^3)*C3 + A*sqrt(y) )/sqrt(mew);

        % Do some checks 
        if(deltaT < TOF)
            psi_low = psi;
        else
            psi_up = psi;
        end

        % Get the new psi value (the trident lookin thing)
        psi = (psi_up+psi_low)/2;

        % Do some more value checking and find new C2 and C3
        if(psi > tall_er_ant)
            C2 = ( 1-cos(sqrt(psi)) )/psi;
            C3 = ( sqrt(psi)-sin(sqrt(psi)) )/(sqrt(psi^3));
        elseif(psi < -tall_er_ant)
            C2 = ( 1-cosh(sqrt(-psi)) )/psi;
            C3 = ( sinh(sqrt(-psi))-sqrt(-psi) )/(sqrt(-psi^3));
        else % probably bad but hey this is what it says to do
            C2 = 1/2;
            C3 = 1/6;
        end
    end
    
    % Get final values, f, g, g_dot
    f = 1 - (y/r_initial);
    g = A*sqrt(y/mew);
    g_dot = 1 - (y/r_final);

    % Now, get the REAL final values
    vel_initial_vec = (pos_vec_final-f*pos_vec_initial)/g;
    vel_final_vec = (pos_vec_final*g_dot-pos_vec_initial)/g;
    return
end