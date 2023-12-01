function [Vol] = CreateVol(deltas)
% CreateVol Computes the volumes of the cells of the domain
% INPUT
% deltas -> spacing between the domain points
% OUTPUT
% Vol -> volumes of the cells (column vector)
Vol = (deltas(2:end-2) + deltas(3:end-1)) / 2;
Vol = [deltas(1) + deltas(2)/2; Vol'; deltas(end) + deltas(end-1)/2];
end
