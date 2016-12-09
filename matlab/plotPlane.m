function plotPlane(point, u, v, xlim, zlim, col)

if nargin < 5
    col = 'g';
end

x = xlim;
z = zlim;
X = [x(1), x(1), x(2), x(2)];
Z = [z(1), z(2), z(2), z(1)];

% Build the matrices to solve for the y coordinates
A = [u(1) v(1); u(3) v(3)];
% Coords of the box of the plane relative to the point
b = [X - point(1); Z - point(3)];
% Factors (a,b) of each vector for the 4 points
factors = A\b;

Y = [u(2) v(2)] * factors + point(2);

h = fill3(X, Z, Y, col);
alpha(h, 0.5);

end
