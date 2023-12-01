function [P] = PhysicsConstants(P)
% PhysicsConstants Sets the physics constants for the parameter structure P
% INPUT
% P -> parameter structure
% OUTPUT
% P -> parameter structure with the physics constants
P.h = 6.62607015e-34;
P.e = 1.602176634e-19;
P.kB = 1.380649e-23;
P.eps0 = 8.854187817e-12;
P.abs0 = 273.15;
P.A = 1.20173e6;
end

