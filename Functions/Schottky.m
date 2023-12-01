function [gamma] = Schottky(E, aT2exp, kBT, beta)
% Schottky Computes the flux emitted from the electrode when a field E is applied
% INPUT
% aT2exp -> 1x2 vector containing coefficients for electrons and holes
% kBT -> Boltzmann constant times the temperature in K
% beta -> coefficient
% E -> 1x2 row vector or matrix of row vectors containing the values of the
% electric field at the interfaces: the first element is E at the left
% electrode (electrons emission), the second element is E at the right
% electrode (holes emission)
% OUTPUT
% gamma -> matrix with the same shape as E in input containing the fluxes (first column
% is referring to electrons emitted from left electrode and second column is referring
% to holes emitted from right electrode)
arguments
    E (:,2) {mustBeNumeric}
    aT2exp (1,2) {mustBeNumeric}
    kBT (1,1) {mustBeNumeric}
    beta (1,1) {mustBeNumeric}
end
gamma = aT2exp .* ( exp( beta*sqrt( abs(E) )/kBT ) - 1 );
end
