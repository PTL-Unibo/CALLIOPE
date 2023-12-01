function [out] = RunExplicit(P, time_instants, options)
% RunExplicit Function to perform a simulation with explicit scheme
% INPUT
% P -> structure containing all the parameters
% time_instants -> vector containing the time instants where the solution
% will be outputted
% options -> structure containing the options for the simulation
% OUTPUT
% out -> structure containing the various results of the simulation

num_instants = length(time_instants);
out.nh = zeros(P.num_points, num_instants);
out.ne = zeros(P.num_points, num_instants);
out.nht = zeros(P.num_points, num_instants);
out.net = zeros(P.num_points, num_instants);
out.rho = zeros(P.num_points, num_instants);
out.phi = zeros(P.num_points+2, num_instants);
out.E = zeros(P.num_points+1, num_instants);
out.J_Sato = zeros(1, num_instants);
out.tout = zeros(1, num_instants);

n = ones(P.num_points, 4) .* P.n_start;
save_flag = true;
save_index = 1;
t = time_instants(1);

start_time_explicit = tic;
while true

    % at this point the number density of all the species at the current time
    % instant are known, we need to calculate the number density of all the
    % species at the next time instant:

    [results] = ExplicitfuncDriftDiffusion(n, P, options);
    
    % Saving the values if the "save_flag" is "true"
    if save_flag
        out.nh(:,save_index) = results.nh;
        out.ne(:,save_index) = results.ne;
        out.nht(:,save_index) = results.nht;
        out.net(:,save_index) = results.net;
        out.rho(:,save_index) = results.rho;
        out.phi(:,save_index) = results.phi;
        out.E(:,save_index) = results.E;
        out.J_Sato(save_index) = results.J_Sato;
        out.tout(save_index) = t; 
        save_index = save_index + 1;
        save_flag = false;
    end

    % Exit the loop if we reached the last time instant
    if save_index > num_instants
        break
    end

    % Computing the time step (dt) that will be used for all the species (minimum of all)
    if options.cfl_type == "advanced"
        if options.source == "On"
            dt = options.cfl / max(max(results.den_for_stab));
        elseif options.source == "Off"
            dt = options.cfl / max(max(results.den_for_stab(:,1:2)));
        end
    elseif options.cfl_type == "simple"
        dt = options.cfl * P.Vol(5) / (2*P.D_h + 1e7*P.mu_h);
    end
    
    % Modifying the dt in order to save at the right time instants
    if t + dt >= time_instants(save_index)
        save_flag = true;
        dt = time_instants(save_index) - t;
    end
    t = t + dt;

    n = n + results.dndt * dt;

    if find(n<0)
        error("Number density became less than 0")
    end

end
out.wct_ppt = toc(start_time_explicit);

end
