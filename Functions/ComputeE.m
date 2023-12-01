function [E] = ComputeE(phi, PhiW, PhiE, deltas)
% ComputeE Function to calculate the electric field
% INPUT
% phi -> column vector or matrix of column vectors containing the electric
% potential
% PhiW -> fixed value of electric potential at the left electrode (1x1)
% PhiE -> fixed value of electric potential at the right electrode (1x1)
% deltas -> spacing between the domain points
% OUTPUT
% E -> column vector or matrix of column vectors containing the electric
% field vaues
arguments
    phi
    PhiW (1,1)
    PhiE (1,1)
    deltas (:,1)
end
E = -(phi(2:end,:) - phi(1:end-1,:)) ./ deltas(2:end-1);
rW = deltas(2) / deltas(1);
rE = deltas(end-1) / deltas(end);
EW = (-phi(1,:) * (rW+1)^2 + (rW+2) * rW * PhiW + phi(2,:)) / (rW * (rW + 1) * deltas(1));
EE = (phi(end,:) * (rE+1)^2 - (rE+2) * rE * PhiE - phi(end-1,:)) / (rE * (rE + 1) * deltas(end));
E = [EW; E; EE];
end
