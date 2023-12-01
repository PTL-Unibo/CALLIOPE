function [x, x_int, x_face] = CreateX(deltas)
% CreateX Creates the x of the domain
% INPUT
% deltas -> spacing between the domain points
% OUTPUT
% x_int -> x coordinates of the points inside the cells
% x -> x coordinates of the points inside the cells plus the coordinates of
% the left and right electrodes
% x_face -> x coordinates of the interfaces between adjacent cells
dim_deltas = length(deltas);
x = zeros(1, dim_deltas+1);
x(1) = 0;
for i = 2:dim_deltas+1
    x(i) = x(i-1) + deltas(i-1);
end
x_int = x(2:end-1);
x_face = (x_int(1:end-1) + x_int(2:end)) / 2;
x_face = [x(1), x_face, x(end)];
end
