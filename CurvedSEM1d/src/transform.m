function [xi, transformed_rho, transformed_mu] = transform(num_spec_el, num_gll, rho_fun, mu_fun)
%Given global points and function handles for the density and rigidity,
%returns the appropriately transformed quantites for optimal solution of
%the wave equation.
syms sym_x
gll_points = double(vpasolve((1 - sym_x^2)*diff(legendreP(num_gll - 1, sym_x)) == 0));
anchor_points = linspace(0, 1, num_spec_el + 1);

global_points_fun = @(x) 2 * pi * (x + 0.5*(1 + gll_points) / num_spec_el);

global_points = arrayfun(global_points_fun, anchor_points(1:end - 1), 'UniformOutput', false);
global_points = uniquetol(cell2mat(global_points), 1 / (num_spec_el * num_gll * 10));

one_over_c = @(x) sqrt(rho_fun(x) ./ mu_fun(x));

int_one_over_c = zeros(size(global_points));
for ii = [1:length(global_points)]
    int_one_over_c(ii) = integral(one_over_c, global_points(1), global_points(ii));
end
F_xi = int_one_over_c(end) ./ global_points(end) .* global_points;
xi = interp1(int_one_over_c, global_points, F_xi);
dxi_dx = int_one_over_c(end) ./ one_over_c(xi) ./ (2*pi);
transformed_rho = dxi_dx .* rho_fun(xi);
transformed_mu = mu_fun(xi) ./ dxi_dx;