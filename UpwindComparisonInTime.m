addpath("Functions\")
clearvars, clc, close all

% Geometry
ns = 16; % number of segments
L = 1; % domain length
ls = L / ns; % segment length
nps = 100; % number of points per segment
np = nps * ns; % total number of points
dx = L / np; % spacing between adjacent points

% Parameters
U = 1;
Amplitude = 9e9;
BaseValue = 1e5;
sigma = 0.75 * ls;

dt_exact = dx / abs(U);
T = L / abs(U);

deltas = ones(1, np+1) * dx;
deltas([1, end]) = deltas([1, end]) / 2;
[x, x_int, x_face] = CreateX(deltas);
delta_x_face = x_face(2:end) - x_face(1:end-1);
Vol = CreateVol(deltas);

n0 = zeros(np,1);
is = @(s) (s-1)*nps+1:s*nps; % function to obtain segment indices
n0(is(3)) = Amplitude * (x_int(is(3)) / ls - 2);
n0(is(4)) = Amplitude * (-x_int(is(4)) / ls + 4);
n0([is(7), is(8)]) = Amplitude;
n0([is(11), is(12), is(13), is(14)]) = Amplitude * exp(-((x_int([is(11), is(12), is(13), is(14)]) - 12*ls)/sigma).^2);
n0 = n0 + BaseValue;

N_exact_matrix = zeros(np); 
for k = 1:np
    N_exact_matrix(:,k) = circshift(n0, sign(U)*(k-1));
end

n_exact = @(t) N_exact_matrix(:, 1 + mod(round(t/dt_exact), np));

time_instants = 0:dt_exact:T;
[tout, nout] = ode45(@(t,n_state)odefunction(n_state,np,U,Vol), time_instants, [n0; n0]);

%%
step = 1;
gif_frame = 15;
create_array = ["gif", "video"];

[plot_options_exact, plot_options_koren, plot_options_upwind] = plot_options_video();

check_folder("media")

for create = create_array

    if create == "video"
        v = VideoWriter(exportpath() + "media\upwind.mp4", "MPEG-4");
        open(v)
    end
    
    fig = figure('WindowState','maximized');
    ax = axes(fig);
    hold on
    set(ax, 'FontSize', 15)
    xlabel('thickness ($\mathrm{m}$)', 'Interpreter','latex')
    ylabel('number density ($\mathrm{m^{-3}}$)', 'Interpreter','latex')
    ax.TickLabelInterpreter = "latex";
    ID = gobjects(3,1);
    ID(1) = plot(x_int(1:step:np), nout(1,1:step:np), plot_options_upwind);
    ID(2) = plot(x_int(1:step:np), nout(1,np+1:step:2*np), plot_options_koren);
    ID(3) = plot(x_int, n_exact(time_instants(1)), plot_options_exact);
    legend('Interpreter','latex', 'Orientation','horizontal')
    grid on
    if create == "video"
        frame = getframe(gcf);
        writeVideo(v,frame)
    elseif create == "gif"
        exportgraphics(gcf,exportpath() + "media\upwind.gif");
    else
        pause(1)
    end
    for k = 2:length(time_instants)
        delete(ID)
        ID(1) = plot(x_int(1:step:np), nout(k,1:step:np), plot_options_upwind);
        ID(2) = plot(x_int(1:step:np), nout(k,np+1:step:2*np), plot_options_koren);
        ID(3) = plot(x_int, n_exact(time_instants(k)), plot_options_exact);
        if create == "video"
            frame = getframe(gcf);
            writeVideo(v,frame)
        elseif create == "gif"
            if mod(k, gif_frame) == 1
                exportgraphics(gcf,exportpath() + "media\upwind.gif",'Append',true);
            end
        else
            pause(0.01)
        end
    end
    
    if create == "video"
        close(v)
    end
    close all

end

%%
leg_x = 0.3;
leg_y = 0.865;
[plot_options_exact, plot_options_koren, plot_options_upwind] = plot_options();

step = 16;

fig = figure;
ax = axes(fig);
hold on
set(ax, 'FontSize', 15)
xlabel('thickness ($\mathrm{m}$)', 'Interpreter','latex')
ylabel('number density ($\mathrm{m^{-3}}$)', 'Interpreter','latex')
ax.TickLabelInterpreter = "latex";
% k = 1;
k = length(tout);
plot(x_int(1:step:end), nout(k,1:step:np), plot_options_upwind)
plot(x_int(1:step:end), nout(k,np+1:step:2*np), plot_options_koren)
plot(x_int, n_exact(time_instants(k)), plot_options_exact)
leg = legend('Interpreter','latex', 'Orientation','horizontal');
set(leg, 'Position', [leg_x, leg_y, ax.Position(1)+ax.Position(3)-leg_x, ax.Position(2)+ax.Position(4)-leg_y])
grid on
box on

%%
check_folder("pdf")
exportgraphics(fig, exportpath() + "pdf\koren_vs_upwind_in_time.pdf")

%%
function [dndt] = odefunction(n, np, U, Vol)
% UPWIND
n_temp_upwind = [n(np); n(1:np); n(1)];
if U > 0
    Gamma_interfaces_upwind = U * n_temp_upwind(1:end-1,:); % if U > 0
elseif U < 0
    Gamma_interfaces_upwind = U * n_temp_upwind(2:end,:); % if U < 0
end

% KOREN
n_temp_koren = [n(end-1); n(end); n(np+1:end); n(np+1); n(np+2)];
a = n_temp_koren(3:end-1) - n_temp_koren(2:end-2);
if U > 0
    b = n_temp_koren(2:end-2) - n_temp_koren(1:end-3); % if U > 0
    koren_computed = KorenMlim(a, b);
    Gamma_interfaces_koren = U * (n_temp_koren(2:end-2,:) + koren_computed); % if U > 0
elseif U < 0
    b = n_temp_koren(4:end) - n_temp_koren(3:end-1); % % if U < 0
    koren_computed = KorenMlim(a, b);
    Gamma_interfaces_koren = U * (n_temp_koren(3:end-1,:) - koren_computed); % if U < 0
end

Gamma_close_upwind = (Gamma_interfaces_upwind(2:end) - Gamma_interfaces_upwind(1:end-1)) ./ Vol;
Gamma_close_koren = (Gamma_interfaces_koren(2:end) - Gamma_interfaces_koren(1:end-1)) ./ Vol;
dndt = - [Gamma_close_upwind; Gamma_close_koren];
end

function [plot_options_exact, plot_options_koren, plot_options_upwind] = plot_options()
    plot_options_exact = struct;
    plot_options_exact.LineStyle = "-";
    plot_options_exact.LineWidth = 1;
    plot_options_exact.Color = [0, 0, 0];
    plot_options_exact.DisplayName = 'analytical';
    
    plot_options_koren = struct;
    plot_options_koren.LineStyle = "-";
    plot_options_koren.LineWidth = 1;
    plot_options_koren.Marker = "o";
    plot_options_koren.MarkerSize = 8;
    plot_options_koren.Color = [1, 0, 1];
    plot_options_koren.DisplayName = 'SOU/KL';
    
    plot_options_upwind = struct;
    plot_options_upwind.LineStyle = "-";
    plot_options_upwind.LineWidth = 1;
    plot_options_upwind.Marker = ".";
    plot_options_upwind.MarkerSize = 15;
    plot_options_upwind.Color = [0, 0, 1];
    plot_options_upwind.DisplayName = 'FOU';
end

function [plot_options_exact, plot_options_koren, plot_options_upwind] = plot_options_video()
    plot_options_exact = struct;
    plot_options_exact.LineStyle = "-";
    plot_options_exact.LineWidth = 1;
    plot_options_exact.Color = [0, 0, 0];
    plot_options_exact.DisplayName = 'analytical';
    
    plot_options_koren = struct;
    plot_options_koren.LineStyle = "-";
    plot_options_koren.LineWidth = 2;
%     plot_options_koren.Marker = "o";
%     plot_options_koren.MarkerSize = 8;
    plot_options_koren.Color = [1, 0, 1];
    plot_options_koren.DisplayName = 'SOU/KL';
    
    plot_options_upwind = struct;
    plot_options_upwind.LineStyle = "-";
    plot_options_upwind.LineWidth = 2;
%     plot_options_upwind.Marker = ".";
%     plot_options_upwind.MarkerSize = 15;
    plot_options_upwind.Color = [0, 0, 1];
    plot_options_upwind.DisplayName = 'FOU';
end
