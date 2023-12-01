function [P] = ParametersFunction(PARAMETER_ID_NAME)
%PARAMETERSFUNCTION contains all the relevant parameter sets used
%   Detailed explanation goes here
arguments
    PARAMETER_ID_NAME char {mustBeMember(PARAMETER_ID_NAME,[ ...
        "LE_ROY",...
        "TEST1"...
        ])}
end

switch PARAMETER_ID_NAME        
    case "LE_ROY"
        % Geometry
        P.L = 4e-4;
        P.num_points = 100;
        P.LW = 0;
        P.LE = 0;
        P.nW = 0;
        P.nE = 0;
        % Material
        P.T = 293.15;
        P.eps_r = 2.2;

        % Essential Parameters
        P.Phi_W = 0;
        P.Phi_E = 4e3;
        P.phih = 1.16;
        P.phie = 1.27;
        P.fix_inj = [0, 0; 0, 0];
        P.n_start = [1e18, 1e18, 1e2, 1e2];
        P.Ndeep = ones(P.num_points,2) .* [6.2e20, 6.2e20];
        
        % Fixed parameters not depending on the electric field 
        P.mu_h = 2e-13;
        P.mu_e = 1e-14;
        P.Bh = 2e-1;
        P.Be = 1e-1;
        P.Dh = 3e-5; %5.8348e-05
        P.De = 3e-5; %1.9133e-04
        P.S0 = 6.4e-22;
        P.S1 = 6.4e-22;
        P.S2 = 6.4e-22;
        P.S3 = 0;

    case "TEST1"
        % Geometry
        P.L = 5e-4;
        P.num_points = 100;
        P.LW = 0;
        P.LE = 0;
        P.nW = 0;
        P.nE = 0;
        % Material
        P.T = 293.15;
        P.eps_r = 2.5;

        % Essential Parameters
        P.Phi_W = 0;
        P.Phi_E = 3e3;
        P.phih = 1.26;
        P.phie = 1.1;
        P.fix_inj = [0, 0; 0, 0];
        P.n_start = [1e19, 1e18, 1e2, 1e2];
        P.Ndeep = ones(P.num_points,2) .* [5e20, 6.2e20];
        
        % Fixed parameters not depending on the electric field 
        P.mu_h = 1e-13;
        P.mu_e = 1e-14;
        P.Bh = 3e-1;
        P.Be = 1e-1;
        P.Dh = 3e-5;
        P.De = 8e-5;
        P.S0 = 6.7e-22;
        P.S1 = 5.4e-22;
        P.S2 = 6.7e-22;
        P.S3 = 0;
end

% Complete P
P = PhysicsConstants(P);
P = DerivedParameters(P);

end
