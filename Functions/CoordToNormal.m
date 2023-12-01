function [norm_x, norm_y] = CoordToNormal(axes_in, x, y)
%CoordToNormal converts a point to normalized coordinates
% INPUT
% axes_in -> handle to axes object
% x -> x coordinate
% y -> y coordinate
% OUTPUT
% norm_x -> normalized x coordinate with respect to figure
% norm_y -> normalized y coordinate with respect to figure
XLIM = axes_in.XLim;
YLIM = axes_in.YLim;
Pos = axes_in.Position;
if axes_in.XScale == "log"
    perc_x = (log(x) - log(XLIM(1))) / (log(XLIM(2)) - log(XLIM(1)));
else
    perc_x = (x - XLIM(1)) / (XLIM(2) - XLIM(1));
end
if axes_in.YScale == "log"
    perc_y = (log(y) - log(YLIM(1))) / (log(YLIM(2)) - log(YLIM(1)));
else
    perc_y = (y - YLIM(1)) / (YLIM(2) - YLIM(1));
end
norm_x = Pos(3) * perc_x + Pos(1); 
norm_y = Pos(4) * perc_y + Pos(2); 
end
