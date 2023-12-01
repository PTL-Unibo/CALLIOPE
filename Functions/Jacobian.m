function [J] = Jacobian(n, Ndeep, B, D, S)
% Jacobian Computes the jacobian of the number density
% n -> 1x4 matrix with the values of the number density for all the
% species
% Ndeep -> 1x2 matrix with the values of the number density of the deep
% traps along the domain
% B -> 1x2 matrix with the trapping coefficients, the first column 
% is referring to holes and the second to electrons
% D -> 1x2 matrix with the detrapping coefficients, the first column 
% is referring to holes and the second to electrons
% S -> 1x4 matrix with the recombination coefficients, the columns are
% referring in order to the recombination between:

arguments
    n (1,4) {mustBeNumeric}
    Ndeep (1,2) {mustBeNumeric}
    B (1,2) {mustBeNumeric}
    D (1,2) {mustBeNumeric}
    S (1,4) {mustBeNumeric}
end

% Initializing the Jacobian matrix
J = zeros(4);

% Naming all the quantities
nh = n(1);
ne = n(2);
nht = n(3);
net = n(4);
Bh = B(1);
Be = B(2);
Dh = D(1);
De = D(2);
S0 = S(1);
S1 = S(2);
S2 = S(3);
S3 = S(4);
Nh = Ndeep(1);
Ne = Ndeep(2);

J(1,1) = -Bh*(1-nht/Nh) - S2*net - S3*ne;
J(1,2) = -S3*nh;
J(1,3) = Bh*(nh/Nh) + Dh;
J(1,4) = -S2*nh;
J(2,1) = -S3*ne;
J(2,2) = -Be*(1-net/Ne) - S1*nht - S3*nh;
J(2,3) = -S1*ne;
J(2,4) = Be*(ne/Ne)+De;
J(3,1) = Bh*(1-nht/Nh);
J(3,2) = -S1*nht;
J(3,3) = -Bh*(nh/Nh) - Dh - S0*net - S1*ne;
J(3,4) = -S0*nht;
J(4,1) = -S2*net;
J(4,2) = Be*(1-net/Ne);
J(4,3) = -S0*net;
J(4,4) = -Be*(ne/Ne) - De - S0*nht - S2*nh;

end
