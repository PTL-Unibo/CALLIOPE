function [EletStat] = EletStat1D(deltas, eps, sparse_flag)
% EletStat1D Creates the elements needed for solving an electrostatic problem in
% a generic 1D domain (the Kelet matrix can be supplied as sparse)
% INPUT
% deltas -> spacing between the domain points
% eps -> dielectric permettivity of the material. It can be a scalar if it is considered 
% constant along the entire domain or a vector if it depends on the position 
% sparse_flag -> if this variable is supplied and is equal to 'sparse' the matrix will be 
% assembled as a sparse matrix, otherwise it will be assembled as a regular matrix
% OUTPUT
% EletStat -> structure containing the elements needed to solve the
% electrostatic problem
arguments
    deltas (1,:) {mustBeNumeric, mustBeReal}
    eps (1,:) {mustBeNumeric, mustBeReal}
    sparse_flag char {mustBeMember(sparse_flag,{'sparse', 'normal'})} = 'normal'
end
num_points = size(deltas, 2) - 1;
EletStat.multRho = (deltas(1:end-1) .* deltas(2:end) ./ (2 * eps))'; 
EletStat.coeffW = deltas(2) / (deltas(2) + deltas(1));
EletStat.coeffE = deltas(end-1) / (deltas(end-1) + deltas(end));
r_minus = deltas(1:end-1)./deltas(2:end);
r_plus = 1 ./ r_minus;
I = [1:num_points, 2:num_points, 1:num_points-1];
J = [1:num_points, 1:num_points-1, 2:num_points];
v = [ones(1,num_points), -1./(1+r_minus(2:end)), -1./(1+r_plus(1:end-1))];
if sparse_flag == "sparse"
    EletStat.Kelet = sparse(I, J, v);
elseif sparse_flag == "normal"
    EletStat.Kelet = zeros(num_points);
    EletStat.Kelet(sub2ind([num_points, num_points],I,J)) = v;
end
end
