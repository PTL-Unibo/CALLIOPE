function [dDdt] = ComputedDdt(E, t, eps)
% ComputedDdt Computes the derivative of vector D over time
% INPUT
% E -> matrix with all the values of electric fiels at the interfaces
% t -> vector containing the time instants at which E is calculated
% eps -> dielectric permettivity
% OUTPUT
% dDdt -> derivative of D over time
arguments
    E {mustBeNumeric}
    t (1,:) {mustBeNumeric}
    eps (1,1) {mustBeNumeric}
end
dEdt = zeros(size(E));
dim_t = size(E,2);
r = zeros(1,dim_t);
delta_t = t(2:end) - t(1:end-1);
r(2:end-1) = delta_t(2:end)./delta_t(1:end-1);
r(1) = delta_t(2)/delta_t(1);
r(end) = delta_t(end-1)/delta_t(end);
dEdt(:,2:end-1) = (E(:,3:end) + (r(2:end-1).^2 - 1).*E(:,2:end-1) - (r(2:end-1).^2).*E(:,1:end-2))./ ...
                (r(2:end-1).*(r(2:end-1) + 1).*delta_t(1:end-1));
dEdt(:,1) = (((r(1)+1)^2)*E(:,2) - r(1)*(r(1)+2)*E(:,1) - E(:,3))/(r(1)*(r(1)+1)*delta_t(1));
dEdt(:,end) = (-((r(end)+1)^2)*E(:,end-1) + r(end)*(r(end)+2)*E(:,end) + E(:,end-2))/(r(end)*(r(end)+1)*delta_t(end));
dDdt = eps * dEdt;
end

