function [out] = Zoom(fig, ax, rect, newax, options)
% Zoom Function to create a zoom in a plot
% INPUT
% fig -> figure handle
% ax -> axes handle
% rect -> struct with fields:
%   rect.cx -> x coordinate of the zoom center
%   rect.cy -> y coordinate of the zoom center
%   rect.w -> width of the zoom box
%   rect.h -> height of the zoom box
% newax -> struct with fields:
%   newax.x -> normalized x coordinate of new axes
%   newax.y -> normalized y coordinate of new axes
%   newax.dx -> normalized dx of new axes
%   newax.dy -> normalized dy of new axes
% start -> string to specify from which vertex the arrow starts
% end -> string to specify to which vertex the arrow arrives
% OUTPUT
% out -> struct with fields:
%   ax -> handle to the new axes created
%   xlim -> vector containing the x limits of the zoom 
%   ylim -> vector containing the y limits of the zoom 
arguments
    fig
    ax
    rect struct
    newax struct
    options.start char {mustBeMember(options.start,{'NW','NE','SW','SE'})} = 'NE'
    options.end char {mustBeMember(options.end,{'NW','NE','SW','SE'})} = 'SW'
end
rectangle(ax, 'Position', [rect.cx-rect.w/2, rect.cy-rect.h/2, rect.w, rect.h])
out.ax = axes(fig, 'Position',[newax.x, newax.y, newax.dx, newax.dy]);
out.xlim = [rect.cx-rect.w/2, rect.cx+rect.w/2];
out.ylim = [rect.cy-rect.h/2, rect.cy+rect.h/2];
switch options.start
    case "NW"
        [x_norm, y_norm] = CoordToNormal(ax, rect.cx-rect.w/2, rect.cy+rect.h/2);
    case "NE"
        [x_norm, y_norm] = CoordToNormal(ax, rect.cx+rect.w/2, rect.cy+rect.h/2);
    case "SW"
        [x_norm, y_norm] = CoordToNormal(ax, rect.cx-rect.w/2, rect.cy-rect.h/2);
    case "SE"
        [x_norm, y_norm] = CoordToNormal(ax, rect.cx+rect.w/2, rect.cy-rect.h/2);
end
switch options.end
    case "NW"
        annotation(fig, 'arrow', [x_norm, newax.x], [y_norm, newax.y + newax.dy])
    case "NE"
        annotation(fig, 'arrow', [x_norm, newax.x + newax.dx], [y_norm, newax.y + newax.dy])
    case "SW"
        annotation(fig, 'arrow', [x_norm, newax.x], [y_norm, newax.y])
    case "SE"
        annotation(fig, 'arrow', [x_norm, newax.x + newax.dx], [y_norm, newax.y])
end
end
