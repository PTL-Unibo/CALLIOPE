function [out] = ExplicitfuncDriftDiffusion(n, P, options)
% ExplicitfuncDriftDiffusion Computes dndt given the current values of n
% INPUT
% n -> matrix containing the number density at the current time step
% P -> structure with all the parameters
% options -> options for the simulation
% OUTPUT
% out.dndt -> time derivative of the number density 
% out.den_for_stab -> npx4 matrix useful to compute the time step dt in an explicit
% solver 

% Solving the Electrostatic problem
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
[Gamma, Gamma_interfaces, den_for_stab_flux] = FluxesExplicit(n(:,1:2), u, P.deltas, P.Vol, P.Diff, BC, options);
[dndt, den_for_stab_omega] = OmegaExplicit(n, P.Ndeep, P.B, P.D, P.S, options);

% Computing the denominator useful for estimating the dt for stability 
den_for_stab = den_for_stab_omega;
den_for_stab(:,1:2) = den_for_stab(:,1:2) + den_for_stab_flux; 

dndt(:,1:2) = dndt(:,1:2) - Gamma;

J_cond = P.e * (Gamma_interfaces(:,1) - Gamma_interfaces(:,2));
J_Sato = -IntegralFunc(J_cond', P.delta_x_face) / P.L;

% Output structure
out.dndt = dndt;
out.den_for_stab = den_for_stab;
out.nh = n(:,1);
out.ne = n(:,2);
out.nht = n(:,3);
out.net = n(:,4);
out.rho = rho;
out.phi = [P.Phi_W; phi; P.Phi_E];
out.E = E;
out.J_Sato = J_Sato;
end
