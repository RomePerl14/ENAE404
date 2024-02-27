% From ENAE301
% TEMPLATE - COPY AND PASTE FOR USE IN FUNCTIONS
% DO NOT CALL THIS FUNCTION YOU GOOBER
function ydot = myodefun(t, y, mew)
    r_mag = norm(y(1:3));
    ydot(1,1) = y(4);
    ydot(2,1) = y(5);
    ydot(3,1) = y(6);
    ydot(4,1) = (-mew/r_mag^3)*y(1);
    ydot(5,1) = (-mew/r_mag^3)*y(2);
    ydot(6,1) = (-mew/r_mag^3)*y(3);
end