function [deltas] = CreateDeltas(LW, LE, nW, nE, Ntot, L)
% CreateDeltas Creates the spacing between the points in the domain
% INPUT
% LW -> length of the thicken spacing region at the left side
% LE -> length of the thicken spacing region at the rigth side
% nW -> number of cells in the left side
% nE -> number of cells in the right side
% Ntot -> total number of cells of the domain
% L -> length of the domain
% OUTPUT
% deltas -> spacing between the domain points
if (LW == 0) || (LE == 0) || (nW == 0) || (nE == 0)
    if (LW == 0) && (LE == 0) && (nW == 0) && (nE == 0)
        d = L / Ntot;
        deltas = ones(1, Ntot+1) * d;
        deltas([1, end]) = deltas([1, end]) / 2;
    else
        error("One or more (but not all) of the input parameters is zero in 'CreateDeltas'")
    end
else
    dW = LW / (nW + 0.5);
    dE = LE / (nE + 0.5);
    deltasW = ones(1, nW+1) * dW;
    deltasW(1) = deltasW(1) / 2;
    deltasE = ones(1, nE+1) * dE;
    deltasE(end) = deltasE(end) / 2;
    nC = Ntot - 2 - nW - nE;
    if nC < 0
        error("The total number of cells is too low in 'CreateDeltas'")
    end
    dC = (L - LW - LE) / (nC + 1);
    deltasC = ones(1, nC+1) * dC;
    deltas = [deltasW, deltasC, deltasE];
end

end
