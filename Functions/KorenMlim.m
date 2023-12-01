function [bphi] = KorenMlim(a, b)
% KorenMlim Function used to compute the Koren flux limiter
% INPUT
% a -> delta in the number density adjacent to the current interface
% b -> delta in the number density adjacent to the upstream interface
% OUTPUT
% bphi -> half the product between b and the flux limiter function

aa = a.*a;  ab = a.*b; % r = a / b

bphi = b; % r > 2.5 (default case)

i = (aa - 2.5*ab) <= 0;
bphi(i) = (b(i) + 2*a(i)) / 6; % r < 2.5

j = (aa - 0.25*ab) <= 0;
bphi(j) = a(j); % r < 0.25 

bphi(ab<=0) = 0; % r < 0 

end
