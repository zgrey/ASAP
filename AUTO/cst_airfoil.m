function [coordU, coordL] = cst_airfoil(l,bl,bu,dt)
%% A function for computing CST airfoil coordinates
% l is a vector discretized upper and lower surface domain from 0 to 1
% bu is a vector of upper surface coefficients in the bernstein polynomial expansion
% bl is a vector of lower surface coefficients in the bernstein polynomial expansion
% dt is the scalar TE offset parameter

% check dimensions of input vectors
if size(l,2) > size(l,1)
    l = l';
end
if size(bu,2) > size(bu,1)
    bu = bu';
end
if size(bl,2) > size(bl,1)
    bl = bl';
end

% "class-function" parameters
N1 = 0.5; N2 = 1;
% compute class function expansion
cf = (l.^N1).*(1-l).^N2;
% respective polynomial degress 
du = length(bu) - 1; dl = length(bl) - 1;
% compute Berstein polynomial basis expansions
Pu = bernsteinMatrix(du,l); Pl = bernsteinMatrix(dl,l);
% compute airfoil coordinates
sU = cf.*(Pu*bu) + dt*l; sL = cf.*(Pl*bl) + dt*l;
% collate into coordinates (function graph)
coordU = [l'; sU']; coordL = [l'; sL'];
