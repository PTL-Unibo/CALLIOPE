function [dndt] = OdefuncDriftDiffusion(t, n_stato, P, options)
% OdefuncDriftDiffusion Computes dndt given the current values of n
% INPUT
% t -> current value of time
% n_stato -> column vector containing the number density
% P -> structure with all the parameters
% options -> options for the simulation
% OUTPUT
% dndt -> first derivative of the number density 
% den_for_stab -> npx4 matrix useful to compute the time step dt in an explicit
% solver 

% Solving the Electrostatic problem
n = reshape(n_stato, [], 4);
rho = sum(n.*[1, -1, 1, -1],2) * P.e;
phi = SolveEletStat(P.EletStat, rho, P.Phi_W, P.Phi_E);
E = ComputeE(phi, P.Phi_W, P.Phi_E, P.deltas);

% Compute mobility
mu = ones(P.num_points+1,2) .* [P.mu_h, P.mu_e];

% The electrodes block the charge
mu([1,end],:) = 0;

% Compute drift velocity 
u = E .* mu .* [1 -1];

% Compute border conditions 
if options.injection == "Schottky"
    gamma = Schottky(E([1, end])', P.aT2exp, P.kBT, P.beta);
    BC = [0, gamma(1); gamma(2), 0];
elseif options.injection == "fixed"
    BC = P.fix_inj;
end

% Compute fluxes and source terms to obtain dndt
[Gamma, ~] = Fluxes(n(:,1:2), u, P.deltas, P.Vol, P.Diff, BC, options);
[dndt] = Omega(n, P.Ndeep, P.B, P.D, P.S, options);

dndt(1:2*P.num_points) = dndt(1:2*P.num_points) - Gamma;
end
