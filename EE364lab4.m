% Code by James Sears 3/5/2024
% Constants
k = 8.9875e9;  % Coulomb's constant

% Define line charge parameters
n_charges = 11;  % Number of charges along the line
length_line = 10;  % Length of the line charge
charges = linspace(-5, 5, n_charges);  % Positions of charges along the line

% Define grid
[x, y, z] = meshgrid(-10:1:10, -10:1:10, -10:1:10);

% Initialize electric field components
Ex = zeros(size(x));
Ey = zeros(size(y));
Ez = zeros(size(z));

% Calculate electric field due to each charge
for i = 1:numel(charges)
    % Distance from each charge to each point in the grid
    dx = x - charges(i);
    dy = y;
    dz = z;

    % Distance squared and cube
    r_squared = dx.^2 + dy.^2 + dz.^2;
    r_cubed = r_squared.^(3/2);

    % Calculate electric field components
    Ex_i = k * dx ./ r_cubed;
    Ey_i = k * dy ./ r_cubed;
    Ez_i = k * dz ./ r_cubed;

    % Superimpose electric field components
    Ex = Ex + Ex_i;
    Ey = Ey + Ey_i;
    Ez = Ez + Ez_i;
end

% Normalize the electric field vectors
E_magnitude = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
Ex = Ex ./ E_magnitude;
Ey = Ey ./ E_magnitude;
Ez = Ez ./ E_magnitude;

% Plot the electric field lines
figure;
quiver3(x, y, z, Ex, Ey, Ez);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Electric Field of a Line Charge');

% Plot for two line charges
figure;
% Define positions of charges
charge_pos1 = [2, 0, 0];  % Position of the first line charge
charge_pos2 = [-2, 0, 0];  % Position of the second line charge
charge1 = 1;  % Charge of the first line (positive)
charge2 = -1;  % Charge of the second line (negative)

% Calculate electric field due to each line charge
% Electric field due to first line charge
dx1 = x - charge_pos1(1);
dy1 = y - charge_pos1(2);
dz1 = z - charge_pos1(3);
r1_squared = dx1.^2 + dy1.^2 + dz1.^2;
r1_cubed = r1_squared.^(3/2);
Ex1 = k * charge1 * dx1 ./ r1_cubed;
Ey1 = k * charge1 * dy1 ./ r1_cubed;
Ez1 = k * charge1 * dz1 ./ r1_cubed;

% Electric field due to second line charge
dx2 = x - charge_pos2(1);
dy2 = y - charge_pos2(2);
dz2 = z - charge_pos2(3);
r2_squared = dx2.^2 + dy2.^2 + dz2.^2;
r2_cubed = r2_squared.^(3/2);
Ex2 = k * charge2 * dx2 ./ r2_cubed;
Ey2 = k * charge2 * dy2 ./ r2_cubed;
Ez2 = k * charge2 * dz2 ./ r2_cubed;

% Superimpose electric field components
Ex_total = Ex1 + Ex2;
Ey_total = Ey1 + Ey2;
Ez_total = Ez1 + Ez2;

% Normalize the electric field vectors
E_total_magnitude = sqrt(Ex_total.^2 + Ey_total.^2 + Ez_total.^2);
Ex_total = Ex_total ./ E_total_magnitude;
Ey_total = Ey_total ./ E_total_magnitude;
Ez_total = Ez_total ./ E_total_magnitude;

% Plot the electric field lines
quiver3(x, y, z, Ex_total, Ey_total, Ez_total);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Electric Field of Two Parallel Line Charges');
