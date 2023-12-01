%%
addpath("Functions\")
clear, clc, close all
np = 100;
U = 5;
D = 0.1;

nL = 5e9;
nR = 2e10;

L = 1;
dx = L / np;
x = ((dx/2):dx:(L-dx/2))';

subx = 1:np;

%%
I = repmat(1:np,1,2);
J = [1:np, 2:np+1];
S = [-ones(1,np), ones(1,np)];
FluxToCell = sparse(I,J,S,np,np+1) / dx;

I = repmat(2:np,1,2);
J = [2:np, 1:np-1];
S = [ones(1,np-1), -ones(1,np-1)];
Kd = -D * sparse(I,J,S,np+1,np) / dx;

Kd(1,1) = -2 * D / dx;
Kd(np+1,np) = +2 * D / dx;

I = 2:np;
if U > 0
    J = 1:np-1;
elseif U < 0
    J = 2:np;
end
S = ones(1,np-1) * U;
Ku = sparse(I,J,S,np+1,np);

Extra = zeros(np+1,1);
Extra(1) = nL * (2*D/dx + U);
Extra(np+1) = -nR * (2*D/dx - U);

A = FluxToCell*(Kd + Ku);
rhs = -FluxToCell * Extra;

K1 = FluxToCell*Kd;
K2 = FluxToCell*Extra;

f = @(n) K1*n + FluxToCell*gammakoren(n, nL, nR, U) + K2;
Jacobian = @(func,x_eval,h_row) (func(repmat(x_eval,size(x_eval')) + diag(h_row)) - func(repmat(x_eval,size(x_eval')))) ./ h_row;

n_upwind = A \ rhs;

n_analytical = (nR - nL*exp(U*L/D) + nL*exp(U*x/D) - nR*exp(U*x/D)) / (1 - exp(U*L/D));

h = 1e-3 * ones(1,np);
n_koren = n_analytical;
r = 1;
while r > 1e-5
    dn = -(Jacobian(f,n_koren,h) \ f(n_koren));
    r = norm(dn);
    n_koren = n_koren + 0.5 * dn;
end

%%
fig1 = figure;
ax = axes(fig1);
hold on
set(ax, 'FontSize', 15)
xlabel('thickness ($\mathrm{m}$)', 'Interpreter','latex')
ylabel('$|$\% relative error$|$', 'Interpreter','latex')
ax.TickLabelInterpreter = "latex";
[~, plot_options_koren, plot_options_upwind] = plot_options();
err_perc_upwind = 100 * abs((n_upwind - n_analytical) ./ n_analytical);
plot(x(subx), err_perc_upwind(subx), plot_options_upwind)
hold on
grid on

err_perc_koren = 100 * abs((n_koren - n_analytical) ./ n_analytical);
plot(x(subx), err_perc_koren(subx), plot_options_koren)

legend('Interpreter','latex', 'Orientation','horizontal', 'Location','northwest')
box on

%%
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\percentage_error.pdf")

%%
[plot_options_exact, plot_options_koren, plot_options_upwind] = plot_options();

fig2 = figure;
ax = axes(fig2);
hold on
set(ax, 'FontSize', 15)
xlabel('thickness ($\mathrm{m}$)', 'Interpreter','latex')
ylabel('number density ($\mathrm{m^{-3}}$)', 'Interpreter','latex')
ax.TickLabelInterpreter = "latex";

plot(ax, x(subx), n_upwind(subx), plot_options_upwind)
plot(ax, x(subx), n_koren(subx), plot_options_koren)
plot(ax, x, n_analytical, plot_options_exact)
grid on
legend('Interpreter','latex', 'Orientation','horizontal')
box on

% Zoom1
rect.cx = 0.985;
rect.cy = 1.2e10;
rect.w = 0.01;
rect.h = 3e9;

newax.x = 0.63;
newax.y = 0.55;
newax.dx = 0.2;
newax.dy = 0.2;

[out1] = Zoom(fig2, ax, rect, newax, "start","NW", "end","NE");
plot(out1.ax, x, n_upwind, plot_options_upwind)
hold on
plot(out1.ax, x, n_koren, plot_options_koren)
plot(out1.ax, x, n_analytical, plot_options_exact)
xlim(out1.xlim)
ylim(out1.ylim)
set(gca, 'FontSize', 12)
grid on

% Zoom2
rect.cx = 0.92;
rect.cy = 0.58e10;
rect.w = 0.08;
rect.h = 1.6e9;

newax.x = 0.23;
newax.y = 0.28;
newax.dx = 0.3;
newax.dy = 0.3;

[out2] = Zoom(fig2, ax, rect, newax, "start","SW", "end","NE");
plot(out2.ax, x, n_upwind, plot_options_upwind)
hold on
plot(out2.ax, x, n_koren, plot_options_koren)
plot(out2.ax, x, n_analytical, plot_options_exact)
xlim(out2.xlim)
ylim(out2.ylim)
set(gca, 'FontSize', 12)
grid on

%%
check_folder("pdf")
exportgraphics(fig2, exportpath() + "pdf\koren_vs_upwind_steadystate.pdf")

%%
function [flux] = gammakoren(n, nL, nR, U)
    nc = size(n,2); 
    n_modified = [nL * ones(1,nc); n; nR * ones(1,nc)];
    a = n_modified(3:end-1,:) - n_modified(2:end-2,:);
    if U > 0
        b = n_modified(2:end-2,:) - n_modified(1:end-3,:); % if U > 0
        koren_computed = KorenMlim(a, b);
        flux = U * (n_modified(2:end-2,:) + koren_computed); % if U > 0
    elseif U < 0
        b = n_modified(4:end,:) - n_modified(3:end-1,:); % % if U < 0
        koren_computed = KorenMlim(a, b);
        flux = U * (n_modified(3:end-1,:) - koren_computed); % if U < 0
    end
    flux = [zeros(1,nc); flux; zeros(1,nc)];
end

function [plot_options_exact, plot_options_koren, plot_options_upwind] = plot_options()
    plot_options_exact = struct;
    plot_options_exact.LineStyle = ":";
    plot_options_exact.LineWidth = 1;
    plot_options_exact.Color = [0, 0, 0];
    plot_options_exact.DisplayName = 'analytical';
    
    plot_options_koren = struct;
    plot_options_koren.LineStyle = "--";
    plot_options_koren.LineWidth = 2;
    plot_options_koren.Color = [1, 0, 1];
    plot_options_koren.DisplayName = 'SOU/KL';
    
    plot_options_upwind = struct;
    plot_options_upwind.LineStyle = "-";
    plot_options_upwind.LineWidth = 2;
    plot_options_upwind.Color = [0, 0, 1];
    plot_options_upwind.DisplayName = 'FOU';
end
