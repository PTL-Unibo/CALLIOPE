function [out] = RunImplicit(P, time_instants, options)
% RunImplicit Function to perform a simulation with ODE
% INPUT
% P -> structure containing all the parameters
% time_instants -> vector containing the time instants where the solution
% will be outputted
% options -> structure containing the options for the simulation
% OUTPUT
% out -> structure containing the various results of the simulation

% Setting initial condition for the number density
n_stato_0 = ones(P.num_points, 4) .* P.n_start;
n_stato_0 = reshape(n_stato_0, [P.num_points*4, 1]);

% Solving with ODE
start_time_ODE = tic;
[out.tout, out.nout] = ode23tb(@(t,n_stato)OdefuncDriftDiffusion(t, n_stato, P, options), time_instants, n_stato_0);
out.wct = toc(start_time_ODE);

% Post processing
start_time_Post_Processing = tic;
if length(out.tout) == length(time_instants)
    [out.nh, out.ne, out.nht, out.net, out.rho, out.phi, out.E, out.J_Sato, out.J_dDdt] = PostProcessing(out.nout, out.tout, P, options);
end
out.ppt = toc(start_time_Post_Processing);

end
