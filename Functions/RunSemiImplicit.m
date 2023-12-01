function [out] = RunSemiImplicit(P, time_instants, options)
% RunSemiImplicit Function to perform a simulation with semi-implicit scheme
% INPUT
% P -> structure containing all the parameters
% time_instants -> vector containing the time instants where the solution
% will be outputted
% options -> structure containing the options for the simulation
% OUTPUT
% out -> structure containing the various results of the simulation

% The Semi-Implicit algorithm is written in Fortran 90
% Here we simply load the output of the Fortran simulation
% (future versions of the code will include a semi-implicit solver in MATLAB)
out.J_Sato = readmatrix('data\data_i.txt');
out.J_Sato = [out.J_Sato(1); out.J_Sato];
out.tout = readmatrix('data\data_t.txt');
end
