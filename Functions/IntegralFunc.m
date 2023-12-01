function [I] = IntegralFunc(f, deltas)
% IntegralFunc integrates f over a domain with points spaced by deltas 
% INPUT
% f -> row vector(1xn) or matrix of row vectors containing the values of the function
% deltas -> vector(n-1) containing the spacing between the points where f is evaluated
% OUTPUT
% I -> column vector containing the value of the integral
arguments
    f {mustBeNumeric}
    deltas (1,:) {mustBeNumeric}
end
I = sum((f(:,1:end-1) + f(:,2:end)) .* deltas,   2) / 2;
end

