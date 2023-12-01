function [out] = Run(P, time_instants, options)
% Run Function to perform a simulation
% INPUT
% P -> structure containing all the parameters
% time_instants -> vector containing the time instants where the solution
% will be outputted
% options -> structure containing options for the simulation
% OUTPUT
% out -> structure containing the various results of the simulation
arguments
    P struct
    time_instants {mustBeNumeric}
    options.time_integration_scheme char {mustBeMember(options.time_integration_scheme,{'explicit', 'semi_implicit', 'implicit'})} = 'implicit'
    options.flux_scheme char {mustBeMember(options.flux_scheme,{'koren', 'upwind'})} = 'upwind'
    options.source char {mustBeMember(options.source, {'On', 'Off'})} = 'On'
    options.injection char {mustBeMember(options.injection, {'Schottky', 'fixed'})} = "Schottky"
    options.cfl {mustBeNumeric} = 1;
    options.cfl_type char {mustBeMember(options.cfl_type, {'simple', 'advanced'})} = "advanced"
end

if options.flux_scheme == "koren"
    if options.time_integration_scheme ~= "implicit"
        warning("Setting the time integration scheme to 'implicit'")
        options.time_integration_scheme = "implicit";
    end
end

switch options.time_integration_scheme
    case "explicit"
        [out] = RunExplicit(P, time_instants, options);
    case "semi_implicit"
        [out] = RunSemiImplicit(P, time_instants, options);
    case "implicit"
        [out] = RunImplicit(P, time_instants, options);
end

end
