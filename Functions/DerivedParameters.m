function [P] = DerivedParameters(P)
% DerivedParameters Updates P with the parameters that depend on the input
% of the simulation but need to be calculated only one time
% INPUT
% P -> parameter structure
% OUTPUT
% P -> parameter structure with derived parameters
P.a = P.A / P.e;
P.deltas = CreateDeltas(P.LW, P.LE, P.nW, P.nE, P.num_points, P.L);
[P.x, P.x_int, P.x_face] = CreateX(P.deltas);
P.delta_x_face = P.x_face(2:end) - P.x_face(1:end-1);
P.Vol = CreateVol(P.deltas);
P.kBT = P.kB * P.T;
P.eps = P.eps_r * P.eps0;
P.beta = sqrt((P.e^3)/(4*pi*P.eps));
P.EletStat = EletStat1D(P.deltas, P.eps, "sparse");
P.Boltz_num = P.e / P.kBT;
P.aT2exp = P.a * (P.T^2) * exp(-[P.phie, P.phih] * P.Boltz_num); 
P.D_h = P.mu_h * P.kBT / P.e;
P.D_e = P.mu_e * P.kBT / P.e;

P.B = ones(P.num_points,2) .* [P.Bh, P.Be];
P.D = ones(P.num_points,2) .* [P.Dh, P.De];
P.S = ones(P.num_points,4) .* [P.S0, P.S1, P.S2, P.S3];
P.Diff = ones(P.num_points+1,2) .* [P.D_h, P.D_e];
end
