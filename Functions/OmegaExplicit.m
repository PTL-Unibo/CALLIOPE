function [omega, den_for_stab] = OmegaExplicit(n, Ndeep, B, D, S, options)
% OmegaExplicit Computes the source terms
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
% OUTPUT
% omega -> column vector containing the source terms for all the species along the domain
% den_for_stab -> npx4 matrix useful to compute the time step dt in an explicit
% solver 

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
den_for_stab = zeros(size(n));

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
    
    % Computing den_for_stab
    den_for_stab(:,1) = Bh.*(1 - nht./Nh) + S2.*net + S3.*ne;
    den_for_stab(:,2) = Be.*(1 - net./Ne) + S1.*nht + S3.*nh;
    den_for_stab(:,3) = Bh.*nh./Nh + Dh + S0.*net + S1.*ne;
    den_for_stab(:,4) = Be.*ne./Ne + De + S0.*nht + S2.*nh;
end
    
end
