function [phi] = SolveEletStat(EletStat, rho, PhiW, PhiE)
% SolveEletStat Solves the electrostatic problem
% INPUT
% EletStat -> structure containing the elements needed to solve the
% electrostatic problem
% rho -> column vector or matrix of column vectors containing the 
% charge density in every cell of the domain
% PhiW -> fixed potential in the left side of the domain
% PhiE -> fixed potential in the right side of the domain
% OUTPUT
% phi -> column vector or matrix of column vectors containing the electric
% potential
rho = rho .* EletStat.multRho;
rho(1,:) = rho(1,:) + EletStat.coeffW * PhiW;
rho(end,:) = rho(end,:) + EletStat.coeffE * PhiE;
phi = EletStat.Kelet\rho;
end
