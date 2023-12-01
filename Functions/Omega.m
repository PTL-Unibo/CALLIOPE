function [omega] = Omega(n, Ndeep, B, D, S, options)
% Omega Computes the source terms
% INPUT
% n -> npx4 matrix with the values of the number density for all the
% species
% Ndeep -> npx2 matrix with the values of the number density of the deep
% traps along the domain
% B -> npx2 matrix with the trapping coefficients, the first column 
% is referring to holes and the second to electrons
% D -> npx2 matrix with the detrapping coefficients, the first column 
% is referring to holes and the second to electrons
% S -> npx4 matrix with the recombination coefficients, the columns are
% referring in order to the recombination between:
% 0) trapped h - trapped e     
% 1) trapped h - mobile e      
% 2) mobile h - trapped e 
% 3) mobile h - mobile e 
% options -> structure containing the options for the simulation
% OUTPUT
% omega -> column vector containing the source terms for all the species along the domain

arguments
    n (:,4) {mustBeNumeric}
    Ndeep (:,2) {mustBeNumeric}
    B (:,2) {mustBeNumeric}
    D (:,2) {mustBeNumeric}
    S (:,4) {mustBeNumeric}
    options struct
end

% Initializing omega and den_for_stab
omega = zeros(size(n));

if options.source == "On"
    % Naming all the quantities
    nh = n(:,1);
    ne = n(:,2);
    nht = n(:,3);
    net = n(:,4);
    Bh = B(:,1);
    Be = B(:,2);
    Dh = D(:,1);
    De = D(:,2);
    S0 = S(:,1);
    S1 = S(:,2);
    S2 = S(:,3);
    S3 = S(:,4);
    Nh = Ndeep(:,1);
    Ne = Ndeep(:,2);
    
    % Assembling the matrix omega following the equations of the model
    Th = Bh.*nh.*(1-nht./Nh);
    Te = Be.*ne.*(1-net./Ne);
    Fh = Dh.*nht;
    Fe = De.*net;
    Rtt = S0.*nht.*net;
    Rtm = S1.*ne.*nht;
    Rmt = S2.*nh.*net;
    Rmm = S3.*nh.*ne;
    omega(:,1) = -Th + Fh - Rmt - Rmm;
    omega(:,2) = -Te + Fe - Rtm - Rmm;
    omega(:,3) = +Th - Fh - Rtt - Rtm;
    omega(:,4) = +Te - Fe - Rtt - Rmt;
end

% Reshape in a single column vector
omega = reshape(omega,[],1);
    
end
