% Fit objective various parameterizations to the NACA airfoil defined below

function phi = coord_obj(x,dv)

% NACA Airfoil
[coordU0,coordL0] = cont_naca4gen(x,0,0,0.12);

% CST fit objective
[coordU, coordL] = CST_airfoil(x',dv(1:length(dv)/2),dv(length(dv)/2+1:end),0);

phi = norm([coordU0(2,:), coordL0(2,:)] - [coordU(2,:), coordL(2,:)]);
