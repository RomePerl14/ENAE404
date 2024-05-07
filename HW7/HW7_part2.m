

%% Q3
mew_sun = 132712*10^6; % km^3/s^2

% Create a 2D matrix containing all of our TOF values
starting_depart_date = ymdhms2jd(2022, 8, 1, 12, 0, 0);
starting_arrival_date = ymdhms2jd(2023, 1, 28, 12, 0, 0);

TOF_matrix(500,500) = 0;
Julian_matrix(500,500) = 0;
for i=0:1:499
    for ii=0:1:499
        TOF_matrix(i+1,ii+1) = ((starting_arrival_date+i) - (starting_depart_date+ii))*24*60*60;
    end
end

test_string = "balls1";

V_infinity_matrix(500,500) = 0;
C3_martix(500,500) = 0;
for i=0:1:499
    for ii=0:1:499
        TOF1 = TOF_matrix(i+1,ii+1);
        [r1_vec1, vel_earth] = findEarth(starting_depart_date+ii);
        [r2_vec1, vel_mars] = findMars(starting_arrival_date+i);
        r1_1 = norm(r1_vec1);
        r2_1 = norm(r2_vec1);
        cos_deltaV1 = (dot(r1_vec1,r2_vec1))/(r1_1*r2_1);
        DM1 = 1; % Given that the thang is short way
        A1 = DM1*sqrt(r1_1*r2_1*(1+cos_deltaV1));
        
        % Weird starter numbers I guess?
        trident1 = 0;
        C2_1 = 1/2;
        C3_1 = 1/6;
        trident_hp1 = 4*(pi^2);
        trident_low1 = -4*(pi^2);
        deltat1 = 0;
        tolerance = 10^-5; % good tolerance
        while(abs(TOF1-deltat1) >= tolerance)
            y1 = r1_1 + r2_1 + (A1*((trident1*C3_1)-1)/sqrt(C2_1));
            if(A1>0 && y1 <0)
                trident_low1 = trident_low1 + (pi/4);
            end
            x1 = sqrt(y1/C2_1);
            deltat1 = (((x1^3)*C3_1)+(A1*sqrt(y1)))/sqrt(mew_sun);
            if deltat1 < TOF1
                trident_low1 = trident1;
            else
                trident_hp1 = trident1;
            end
            trident1 = (trident_hp1 + trident_low1)/2;
            if trident1 > tolerance
                C2_1 = (1-cos(sqrt(trident1)))/trident1;
                C3_1 = (sqrt(trident1)-sin(sqrt(trident1)))/sqrt(trident1^3);
            elseif trident1 < -tolerance
                C2_1 = (1-cosh(sqrt(-trident1)))/trident1;
                C3_1 = (sinh(sqrt(-trident1))-sqrt(-trident1))/sqrt((-trident1)^3);
            else
                C2_1 = 1/2;
                C3_1 = 1/6;
            end
        end
        % Now find other stuff (what is going ON with these variables maine)
        f1 = 1 - (y1/r1_1);
        g1 = A1*sqrt(y1/mew_sun);
        g_dot1 = 1 - (y1/r2_1);
        v1_vec1 = (r2_vec1 - (f1*r1_vec1))/g1; % Initial Vel
        v2_vec1 = ((g_dot1*r2_vec1) - r1_vec1)/g1; % Final vel
        % Fist case, short way
        V_infinity_earth1 = v1_vec1 - vel_earth;
        V_infinity_mars1 = v2_vec1 - vel_mars;
        % Get the magnitude:
        val1 = norm(V_infinity_earth1) + norm(V_infinity_mars1);

        % -- Start again
        TOF1 = TOF_matrix(i+1,ii+1);
        [r1_vec1, vel_earth] = findEarth(starting_depart_date+ii);
        [r2_vec1, vel_mars] = findMars(starting_arrival_date+i);
        r1_1 = norm(r1_vec1);
        r2_1 = norm(r2_vec1);
        cos_deltaV1 = (dot(r1_vec1,r2_vec1))/(r1_1*r2_1);
        DM1 = -1; % Given that the thang is short way
        A1 = DM1*sqrt(r1_1*r2_1*(1+cos_deltaV1));
        
        % Weird starter numbers I guess?
        trident1 = 0;
        C2_1 = 1/2;
        C3_1 = 1/6;
        trident_hp1 = 4*(pi^2);
        trident_low1 = -4*(pi^2);
        deltat1 = 0;
        tolerance = 10^-6; % good tolerance
        while(abs(TOF1-deltat1) >= tolerance)
            y1 = r1_1 + r2_1 + (A1*(trident1*C3_1-1)/sqrt(C2_1));
            if(A1>0 && y1 <0)
                trident_low1 = trident_low1 + (pi/4);
            end
            x1 = sqrt(y1/C2_1);
            deltat1 = ((x1^3)*C3_1+A1*sqrt(y1))/sqrt(mew_sun);
            if deltat1 < TOF1
                trident_low1 = trident1;
            else
                trident_hp1 = trident1;
            end
            trident1 = (trident_hp1 + trident_low1)/2;
            if trident1 > tolerance
                C2_1 = (1-cos(sqrt(trident1)))/trident1;
                C3_1 = (sqrt(trident1)-sin(sqrt(trident1)))/sqrt(trident1^3);
            elseif trident1 < -tolerance
                C2_1 = (1-cosh(sqrt(-trident1)))/trident1;
                C3_1 = (sinh(sqrt(-trident1))-sqrt(-trident1))/sqrt((-trident1)^3);
            else
                C2_1 = 1/2;
                C3_1 = 1/6;
            end 
            fprintf("im stuck + " + int2str(abs(TOF1-deltat1)) + "\n")
        end
        
        % Now find other stuff (what is going ON with these variables maine)
        f1 = 1 - (y1/r1_1);
        g1 = A1*sqrt(y1/mew_sun);
        g_dot1 = 1 - (y1/r2_1);
        v1_vec1 = (r2_vec1 - (f1*r1_vec1))/g1; % Initial Vel
        v2_vec1 = ((g_dot1*r2_vec1) - r1_vec1)/g1; % Final vel
        % Second case, long way
        V_infinity_earth2 = v1_vec2 - vel_earth2;
        V_infinity_mars2 = v2_vec2 - vel_mars2;
        % Get the magnitude:
        val2 = norm(V_infinity_earth2) + norm(V_infinity_mars2);

        % now find which mag sum is smaller, and keep it:
        if (val1 < val2)
            C3_REAL= (norm(V_infinity_earth1))^2;
            V_infinity_REAL = norm(V_infinity_mars1);
        elseif (val2 < val1)
            C3_REAL = (norm(V_infinity_earth2))^2;
            V_infinity_REAL = norm(V_infinity_mars2);
        else % for some reason if either conditional doesn't occur, just keep the short way
            C3_REAL= (norm(V_infinity_earth1))^2;
            V_infinity_REAL = norm(V_infinity_mars1);
        end
        V_infinity_matrix(i+1,ii+1) = V_infinity_REAL;
        C3_martix(i+1,ii+1) = C3_REAL;

    end

end

hold on
contour(TOF_matrix, "-k", DisplayName="TOF")
contour(V_infinity_matrix, "-r", DisplayName="V_i_n_f_i_n_i_t_y")
contour(C3_martix, "-b", DisplayName="C3")


