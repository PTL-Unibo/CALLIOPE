function [J_cond] = ComputeJCond(nh, ne, E, Diff_h, Diff_e, u_h, u_e, deltas, aT2exp, kBT, beta, e, options)
% ComputeJCond Computes the current density due to conduction
% INPUT
% nh -> holes number density 
% ne -> electrons number density 
% E -> matrix with the electric field at all the interfaces 
% Diff_h -> matrix with the diffusion coefficients for holes
% Diff_e -> matrix with the diffusion coefficients for electrons 
% u_h -> velocity for holes 
% u_e -> velocity for electrons 
% deltas -> spacing between the domain points 
% aT2exp -> 1x2 vector with the coefficients for Schottky emission 
% kBT -> Boltzmann constant times temperature
% beta -> constant 
% e -> electron elementary charge
% options -> options for the simulation
% OUTPUT
% J_cond -> conduction current density at all the interfaces for all time
% instants

J_cond_h = zeros(size(E));
J_cond_e = zeros(size(E));

if options.flux_scheme == "upwind" 
    J_diffusion_h = -Diff_h(2:end-1, :) .* (nh(2:end,:) - nh(1:end-1,:)) ./ deltas(2:end-1)';
    J_diffusion_e = -Diff_e(2:end-1, :) .* (ne(2:end,:) - ne(1:end-1,:)) ./ deltas(2:end-1)';
    Umax_h = max(0, u_h);
    umin_h = min(0, u_h);
    Umax_e = max(0, u_e);
    umin_e = min(0, u_e);
    J_drift_h = nh(1:end-1,:).*Umax_h(2:end-1,:) + nh(2:end,:).*umin_h(2:end-1,:);
    J_drift_e = ne(1:end-1,:).*Umax_e(2:end-1,:) + ne(2:end,:).*umin_e(2:end-1,:);
    J_cond_h(2:end-1,:) = J_diffusion_h + J_drift_h;
    J_cond_e(2:end-1,:) = J_diffusion_e + J_drift_e;
elseif options.flux_scheme == "koren"
    J_diffusion_h = -Diff_h(2:end-1, :) .* (nh(2:end,:) - nh(1:end-1,:)) ./ deltas(2:end-1)';
    J_diffusion_e = -Diff_e(2:end-1, :) .* (ne(2:end,:) - ne(1:end-1,:)) ./ deltas(2:end-1)';
    Upos_h = u_h >= 0;
    uneg_h = u_h < 0;
    Umax_h = u_h.*Upos_h;
    umin_h = u_h.*uneg_h;
    Upos_e = u_e >= 0;
    uneg_e = u_e < 0;
    Umax_e = u_e.*Upos_e;
    umin_e = u_e.*uneg_e;

    nh_modified = [zeros(1,size(nh,2)); nh; zeros(1,size(nh,2))];
    DeltaN_h = nh_modified(2:end,:) - nh_modified(1:end-1,:);
    DeltaN_h_to_koren = Upos_h(2:end-1,:).*DeltaN_h(1:end-2,:) + uneg_h(2:end-1,:).*DeltaN_h(3:end,:);
    koren_computed_h = KorenMlim(DeltaN_h(2:end-1,:), DeltaN_h_to_koren);
    J_drift_h = Umax_h(2:end-1,:).*(nh_modified(2:end-2,:) + koren_computed_h) +...
                umin_h(2:end-1,:).*(nh_modified(3:end-1,:) - koren_computed_h);

    ne_modified = [zeros(1,size(ne,2)); ne; zeros(1,size(ne,2))];
    DeltaN_e = ne_modified(2:end,:) - ne_modified(1:end-1,:);
    DeltaN_e_to_koren = Upos_e(2:end-1,:).*DeltaN_e(1:end-2,:) + uneg_e(2:end-1,:).*DeltaN_e(3:end,:);
    koren_computed_e = KorenMlim(DeltaN_e(2:end-1,:), DeltaN_e_to_koren);
    J_drift_e = Umax_e(2:end-1,:).*(ne_modified(2:end-2,:) + koren_computed_e) +...
                umin_e(2:end-1,:).*(ne_modified(3:end-1,:) - koren_computed_e);

    J_cond_h(2:end-1,:) = J_diffusion_h + J_drift_h;
    J_cond_e(2:end-1,:) = J_diffusion_e + J_drift_e;
end

E_to_schottky = [E(1,:)', E(end,:)'];
J_schottky = Schottky(E_to_schottky, aT2exp, kBT, beta);
J_cond_e(1,:) = J_schottky(:,1);
J_cond_h(end,:) = -J_schottky(:,2);

J_cond = e * (J_cond_h - J_cond_e);
end
