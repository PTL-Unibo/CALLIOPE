%% SATO
addpath('Functions')
clearvars, clc, close all

P = ParametersFunction("LE_ROY");

time_instants = [0, logspace(0,5, 60)];

[out] = Run(P, time_instants);

%% SATO GRAPH

plot_options_sato.LineWidth = 2;
plot_options_sato.DisplayName = 'Sato';
plot_options_sato.LineStyle = '-';

plot_options_JdDdt.LineWidth = 2;
plot_options_JdDdt.DisplayName = '$J + \frac{\partial D}{\partial t}$';
plot_options_JdDdt.LineStyle = 'none';
plot_options_JdDdt.Marker = '.';
plot_options_JdDdt.Markersize = 20;

fig1 = figure();
ax1 = axes(fig1);
loglog(ax1, out.tout, out.J_Sato, plot_options_sato)
hold on
loglog(ax1, out.tout, out.J_dDdt, plot_options_JdDdt)

% title('Sato vs. $J + \frac{\partial D}{\partial t}$', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
polarization_current_plot()

rect.cx = 45816; rect.cy = 1.3e-10; rect.w = 4e4; rect.h = 4e-11;
newax.x = 0.52; newax.y = 0.42; newax.dx = 0.3; newax.dy = 0.3;
[out_zoom] = Zoom(fig1, ax1, rect, newax, 'start','NE', 'end','SE');

loglog(out_zoom.ax, out.tout, out.J_Sato, plot_options_sato)
hold on
loglog(out_zoom.ax, out.tout, out.J_dDdt, plot_options_JdDdt)
grid on
set(gca, 'FontSize', 12)

out_zoom.ax.XLim = out_zoom.xlim;
out_zoom.ax.YLim = out_zoom.ylim;

%% SATO EXPORT 
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\sato.pdf")

%% COMPUTATIONAL TIME
addpath('Functions')
clearvars, clc, close all
P = ParametersFunction("LE_ROY");

n = 20;
max_time_array = logspace(1, 5, n);
t_out = zeros(n,2);

i = 1;
for t_max = max_time_array
    time_instants = logspace(0, log10(t_max), 100);
    [out_implicit] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit");
    comput_time_implicit = out_implicit.wct + out_implicit.ppt;
    %fprintf("Implicit: %f\n", comput_time_implicit)
    [out_explicit] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",0.8);
    comput_time_explicit = out_explicit.wct_ppt;
    %fprintf("Explicit: %f\n", comput_time_explicit)
    t_out(i,:) = [comput_time_implicit, comput_time_explicit];
    i = i + 1;
end

%% COMPUTATIONAL TIME GRAPH 

plot_options_comp_exp.LineWidth = 2;
plot_options_comp_exp.DisplayName = 'Explicit';
plot_options_comp_exp.LineStyle = ':';
plot_options_comp_exp.Color = [150, 100, 0]/255;

plot_options_comp_imp.LineWidth = 2;
plot_options_comp_imp.DisplayName = 'Implicit';
plot_options_comp_imp.LineStyle = '-';
plot_options_comp_imp.Color = [1, 0, 0];

fig1 = figure();
ax1 = axes(fig1);
semilogx(ax1, max_time_array, t_out(:,1), plot_options_comp_imp)
hold on
semilogx(max_time_array, t_out(:,2), plot_options_comp_exp)
%title('Computational time', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
computational_time_plot()

%% COMPUTATIONAL TIME EXPORT 
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\computational_time.pdf")

%% KOREN VS UPWIND
addpath('Functions')
clearvars, clc, close all
P = ParametersFunction("LE_ROY");

time_instants = [0, logspace(0,5, 99)];

[out_upwind] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit");
fprintf("Upwind: %f\n", out_upwind.wct + out_upwind.ppt)
time_func_upwind = @() Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit");
[out_koren] = Run(P, time_instants, "flux_scheme","koren", "time_integration_scheme","implicit");
fprintf("Koren: %f\n", out_koren.wct + out_koren.ppt)
time_func_koren = @() Run(P, time_instants, "flux_scheme","koren", "time_integration_scheme","implicit");

% fprintf("timeit upwind: \n") % 1.3028 s
% timeit(time_func_upwind)
% fprintf("timeit koren: \n") % 11.0638 s
% timeit(time_func_koren)

%% KOREN VS UPWIND GRAPH 
% Point of maximum percentage error = 10.964911 %: (7543.120063, 2.795912e-10)

[plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren();

fig1 = figure();
ax1 = axes(fig1);
loglog(ax1, out_upwind.tout, out_upwind.J_Sato, plot_options_upwind);
hold on
loglog(ax1, out_koren.tout, out_koren.J_Sato, plot_options_koren);

[max_err_perc, i] = max(abs(out_upwind.J_Sato - out_koren.J_Sato)*100 ./ out_upwind.J_Sato);
t_max_err_perc = out_upwind.tout(i);
J_max_err_perc_up = out_upwind.J_Sato(i);
J_max_err_perc_kor = out_koren.J_Sato(i);
J_max_err_perc = (J_max_err_perc_up + J_max_err_perc_kor) / 2;
fprintf('Point of maximum percentage error = %f %%: (%f, %e)\n', max_err_perc, t_max_err_perc, J_max_err_perc)

%title('Koren vs. Upwind', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
polarization_current_plot();

rect.cx = t_max_err_perc; rect.cy = J_max_err_perc; rect.w = 2000; rect.h = 5e-11; 
newax.x = 0.5; newax.y = 0.45; newax.dx = 0.3; newax.dy = 0.3;
out_zoom = Zoom(fig1, ax1, rect, newax, "start","SE", "end","SE");

loglog(out_zoom.ax, out_upwind.tout, out_upwind.J_Sato, plot_options_upwind);
hold on
loglog(out_zoom.ax, out_koren.tout, out_koren.J_Sato, plot_options_koren);
plot(out_zoom.ax, t_max_err_perc, J_max_err_perc_up, '.', "Color",plot_options_upwind.Color, 'MarkerSize',20)
plot(out_zoom.ax, t_max_err_perc, J_max_err_perc_kor, '.', "Color",plot_options_koren.Color, 'MarkerSize',20)
grid on
set(gca, 'FontSize', 12)
out_zoom.ax.XTick = [7000, 7500, 8000];

out_zoom.ax.XLim = out_zoom.xlim;
out_zoom.ax.YLim = out_zoom.ylim;

%% KOREN VS UPWIND EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_vs_upwind.pdf")

%% CHARGE DENSITY GRAPH
% Point of maximum error = 0.493720: (0.000394, 6.060046e+00)

[plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren();

% [max_col, ~] = max(abs(out_upwind.rho - out_koren.rho)./ out_upwind.rho);
% [~, err_max_time_index] = max(max_col);
err_max_time_index = 100;

fig1 = figure();
ax1 = axes(fig1);
plot(ax1, P.x_int, out_upwind.rho(:,err_max_time_index), plot_options_upwind);
hold on
plot(ax1, P.x_int, out_koren.rho(:,err_max_time_index), plot_options_koren);

%title('Charge Density', Interpreter="latex")
ax1.TickLabelInterpreter = "latex";
charge_density_plot()

rho_upwind = out_upwind.rho(:,err_max_time_index);
rho_koren = out_koren.rho(:,err_max_time_index);
[max_err, i] = max(abs(rho_upwind - rho_koren));
x_max_err = P.x_int(i);
rho_max_err_up = out_upwind.rho(i,err_max_time_index);
rho_max_err_kor = out_koren.rho(i,err_max_time_index);
rho_max_err = (rho_max_err_up + rho_max_err_kor) / 2;
fprintf('Point of maximum error = %f: (%f, %e)\n', max_err, x_max_err, rho_max_err)
rect.cx = x_max_err; rect.cy = rho_max_err; rect.w = 5e-6; rect.h = 1.2; 

% rect.cx = 3.4e-4; rect.cy = 0.9; rect.w = 2e-5; rect.h = 1.2; 

newax.x = 0.48; newax.y = 0.62; newax.dx = 0.24; newax.dy = 0.28;
out_zoom = Zoom(fig1, ax1, rect, newax, "start","NE", "end","NE");

plot(out_zoom.ax, P.x_int, out_upwind.rho(:,err_max_time_index), plot_options_upwind);
hold on
plot(out_zoom.ax, P.x_int, out_koren.rho(:,err_max_time_index), plot_options_koren);
plot(out_zoom.ax, x_max_err, rho_max_err_up, '.', "Color",plot_options_upwind.Color, 'MarkerSize',20)
plot(out_zoom.ax, x_max_err, rho_max_err_kor, '.', "Color",plot_options_koren.Color, 'MarkerSize',20)
grid on
set(gca, 'FontSize', 12)

out_zoom.ax.XLim = out_zoom.xlim;
out_zoom.ax.YLim = out_zoom.ylim;

%% CHARGE DENSITY EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_vs_upwind_rho.pdf")

%% ELECTRONS NUMBER DENSITY GRAPH

[plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren();
[plot_options_upwind, plot_options_koren] = ...
    extra_plot_options_upwind_koren(plot_options_upwind, plot_options_koren);


% [max_col, ~] = max(abs(out_upwind.ne - out_koren.ne)./ out_upwind.ne);
% [~, err_max_time_index] = max(max_col);
err_max_time_index = 100;

fig1 = figure();
ax1 = axes(fig1);
semilogy(ax1, P.x_int, out_upwind.ne(:,err_max_time_index), plot_options_upwind);
hold on
semilogy(ax1, P.x_int, out_koren.ne(:,err_max_time_index), plot_options_koren);

%title('Electrons', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
number_density_plot()

%% ELECTRONS NUMBER DENSITY EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_vs_upwind_ne.pdf")

%% HOLES NUMBER DENSITY GRAPH

[plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren();
[plot_options_upwind, plot_options_koren] = ...
    extra_plot_options_upwind_koren(plot_options_upwind, plot_options_koren);

% [max_col, ~] = max(abs(out_upwind.nh - out_koren.nh)./ out_upwind.nh);
% [~, err_max_time_index] = max(max_col);
err_max_time_index = 100;

fig1 = figure();
ax1 = axes(fig1);
semilogy(ax1, P.x_int, out_upwind.nh(:,err_max_time_index), plot_options_upwind);
hold on
semilogy(ax1, P.x_int, out_koren.nh(:,err_max_time_index), plot_options_koren);

%title('Holes', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
number_density_plot()

%% HOLES NUMBER DENSITY EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_vs_upwind_nh.pdf")

%% TRAPPED ELECTRONS NUMBER DENSITY GRAPH

[plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren();
[plot_options_upwind, plot_options_koren] = ...
    extra_plot_options_upwind_koren(plot_options_upwind, plot_options_koren);

% [max_col, ~] = max(abs(out_upwind.net - out_koren.net)./ out_upwind.net);
% [~, err_max_time_index] = max(max_col);
err_max_time_index = 100;

fig1 = figure();
ax1 = axes(fig1);
semilogy(ax1, P.x_int, out_upwind.net(:,err_max_time_index), plot_options_upwind);
hold on
semilogy(ax1, P.x_int, out_koren.net(:,err_max_time_index), plot_options_koren);

%title('Trapped Electrons', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
number_density_plot()

%% TRAPPED ELECTRONS NUMBER DENSITY EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_vs_upwind_net.pdf")

%% TRAPPED HOLES NUMBER DENSITY GRAPH

[plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren();
[plot_options_upwind, plot_options_koren] = ...
    extra_plot_options_upwind_koren(plot_options_upwind, plot_options_koren);

% [max_col, ~] = max(abs(out_upwind.nht - out_koren.nht)./ out_upwind.nht);
% [~, err_max_time_index] = max(max_col);
err_max_time_index = 100;

fig1 = figure();
ax1 = axes(fig1);
semilogy(ax1, P.x_int, out_upwind.nht(:,err_max_time_index), plot_options_upwind);
hold on
semilogy(ax1, P.x_int, out_koren.nht(:,err_max_time_index), plot_options_koren);

%title('Trapped Holes', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
number_density_plot()

%% TRAPPED HOLES NUMBER DENSITY EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_vs_upwind_nht.pdf")

%% IMPLICIT VS EXPLICIT
addpath('Functions')
clearvars, clc, close all
P = ParametersFunction("LE_ROY");

time_instants = [0, logspace(0,5, 99)];

[out_implicit] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit");
time_func_implicit = @() Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit");
fprintf("Implicit: %f\n", out_implicit.wct + out_implicit.ppt)
[out_explicit] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",0.8);
time_func_explicit = @() Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",0.8);
fprintf("Explicit: %f\n", out_explicit.wct_ppt)

% fprintf("timeit implicit: \n") % 1.3069
% timeit(time_func_implicit)
% fprintf("timeit explicit: \n") % 9.1859
% timeit(time_func_explicit)

%% IMPLICIT VS EXPLICIT GRAPH
% Point of maximum percentage error = 23.937940 %: (21.209509, 5.993979e-09)

plot_options_implicit.LineWidth = 2;
plot_options_implicit.Color = [1, 0, 0];
plot_options_implicit.DisplayName = 'Implicit';

plot_options_explicit.LineWidth = 2;
plot_options_explicit.Color = [150, 100, 0]/255;
plot_options_explicit.DisplayName = 'Explicit';
plot_options_explicit.LineStyle = ':';

fig1 = figure();
ax1 = axes(fig1);
loglog(ax1, out_implicit.tout, out_implicit.J_Sato, plot_options_implicit);
hold on
loglog(ax1, out_explicit.tout, out_explicit.J_Sato, plot_options_explicit);

[max_err_perc, i] = max(abs(out_implicit.J_Sato - out_explicit.J_Sato')*100 ./ out_implicit.J_Sato);
t_max_err_perc = out_implicit.tout(i);
J_max_err_perc_imp = out_implicit.J_Sato(i);
J_max_err_perc_exp = out_explicit.J_Sato(i);
J_max_err_perc = (J_max_err_perc_exp + J_max_err_perc_imp) / 2;
fprintf('Point of maximum percentage error = %f %%: (%f, %e)\n', max_err_perc, t_max_err_perc, J_max_err_perc)

%title('Implicit vs. Explicit', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
polarization_current_plot();

rect.cx = t_max_err_perc; rect.cy = J_max_err_perc; rect.w = 10; rect.h = 3e-9; 
newax.x = 0.5; newax.y = 0.45; newax.dx = 0.3; newax.dy = 0.3;
out_zoom = Zoom(fig1, ax1, rect, newax, "start","SE", "end","SW");

loglog(out_zoom.ax, out_implicit.tout, out_implicit.J_Sato, plot_options_implicit);
hold on
loglog(out_zoom.ax, out_explicit.tout, out_explicit.J_Sato, plot_options_explicit);
plot(out_zoom.ax, t_max_err_perc, J_max_err_perc_imp, '.', "Color",plot_options_implicit.Color, 'MarkerSize',20)
plot(out_zoom.ax, t_max_err_perc, J_max_err_perc_exp, '.', "Color",plot_options_explicit.Color, 'MarkerSize',20)
grid on
set(gca, 'FontSize', 12)

out_zoom.ax.XLim = out_zoom.xlim;
out_zoom.ax.YLim = out_zoom.ylim;

%% IMPLICIT VS EXPLICIT EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\implicit_vs_explicit.pdf")

%% IMPLICIT VS SEMI-IMPLICIT
addpath('Functions')
clearvars, clc, close all
P = ParametersFunction("LE_ROY");

[out_semiimplicit] = Run(struct, zeros(1,10), "time_integration_scheme","semi_implicit");

[out_implicit] = Run(P, out_semiimplicit.tout, "flux_scheme","upwind", "time_integration_scheme","implicit");

%% IMPLICIT VS SEMI-IMPLICIT GRAPH 
% Point of maximum percentage error = 23.544915 %: (7355.508000, 2.652013e-10)

plot_options_implicit.LineWidth = 2;
plot_options_implicit.Color = [1, 0, 0];
plot_options_implicit.DisplayName = 'Implicit';

plot_options_semiimplicit.LineWidth = 2;
plot_options_semiimplicit.Color = [255, 69, 0]/255;
plot_options_semiimplicit.DisplayName = 'Semi-Implicit';
plot_options_semiimplicit.LineStyle = '--';

[max_err_perc, i] = max(abs(out_implicit.J_Sato - out_semiimplicit.J_Sato)*100 ./ out_implicit.J_Sato);
t_max_err_perc = out_semiimplicit.tout(i);
J_max_err_perc_imp = out_implicit.J_Sato(i);
J_max_err_perc_semiimp = out_semiimplicit.J_Sato(i);
J_max_err_perc = (J_max_err_perc_semiimp + J_max_err_perc_imp) / 2;
fprintf('Point of maximum percentage error = %f %%: (%f, %e)\n', max_err_perc, t_max_err_perc, J_max_err_perc)

fig1 = figure();
ax1 = axes(fig1);
imp = loglog(ax1, out_implicit.tout, out_implicit.J_Sato, plot_options_implicit);
hold on
semi = loglog(ax1, out_semiimplicit.tout, out_semiimplicit.J_Sato, plot_options_semiimplicit);

%title('Implicit vs. Semi-Implicit', 'Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
polarization_current_plot();

rect.cx = t_max_err_perc; rect.cy = J_max_err_perc; rect.w = 800; rect.h = 1.2e-10; 
newax.x = 0.5; newax.y = 0.42; newax.dx = 0.3; newax.dy = 0.3;
out_zoom = Zoom(fig1, ax1, rect, newax, "start","SE", "end","SE");

loglog(out_zoom.ax, out_implicit.tout, out_implicit.J_Sato, plot_options_implicit);
hold on
loglog(out_zoom.ax, out_semiimplicit.tout, out_semiimplicit.J_Sato, plot_options_semiimplicit);
plot(out_zoom.ax, t_max_err_perc, J_max_err_perc_imp, '.', "Color",plot_options_implicit.Color, 'MarkerSize',20)
plot(out_zoom.ax, t_max_err_perc, J_max_err_perc_semiimp, '.', "Color",plot_options_semiimplicit.Color, 'MarkerSize',20)
grid on
set(gca, 'FontSize', 12)

out_zoom.ax.XLim = out_zoom.xlim;
out_zoom.ax.YLim = out_zoom.ylim;

%% IMPLICIT VS SEMI-IMPLICIT EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\implicit_vs_semiimplicit.pdf")

%% DOMAIN GRAPH
clearvars, clc, close all
L = 0;
U = 5;
xp = [1, 2, 3, 4];
offset_x = -0.15;
offset_y = 0.42;
offset_x2 = -0.18;
offset_y2 = 1.3;

fig1 = figure();
ax1 = axes(fig1);
hold on
for i = 1:length(xp)
    plot(xp(i), 0, 'k.', 'MarkerSize',25)
    eval("text(xp(i)+offset_x, offset_y, '$n_{i " + my_num2str(i-2) +"}$', 'Interpreter','latex', 'FontSize',15);")
end
for i = 1:length(xp)
    plot([xp(i), xp(i)]-0.5, [-1,1], 'k--', 'LineWidth',0.5)
end
plot([xp(i), xp(i)]+0.5, [-1,1], 'k--', 'LineWidth',0.5)
plot([L,U],[0,0],'k','LineWidth',0.2)
xlim([L, U])
ylim([-1, 1]*4)

text((L + U)/2 + offset_x2, offset_y2, '$i + \frac{1}{2}$', 'Interpreter','latex', 'FontSize',12);

set(gca, 'Visible', 'off')

%% DOMAIN EXPORT
check_folder("pdf")
exportgraphics(ax1, exportpath() + "pdf\domain.pdf")

%% KOREN FUNCTION GRAPH
clearvars, clc, close all
low_lim = -0.5;
up_lim = 3;
lim1 = [-0.5, 0];
lim2 = [-0, 0.25];
lim3 = [0.25, 2.5];
lim4 = [2.5, 3];
num = 200;

col_array = ["#7E2F8E", "#0072BD", "#D95319", "#7E2F8E"];
style_array = ["-", "-", "-", "-"];

r1 = linspace(lim1(1), lim1(2), num);
r2 = linspace(lim2(1), lim2(2), num);
r3 = linspace(lim3(1), lim3(2), num);
r4 = linspace(lim4(1), lim4(2), num);

phi = @(r) max(0, min(2*r, min(2, (1+2*r)/3)));

fig1 = figure();
ax1 = axes(fig1);
hold on
for i = 1:4
    eval("plot(ax1, r" + num2str(i) + ", phi(r" + num2str(i) + "), 'LineWidth',3," + ...
        "'LineStyle', style_array(" + num2str(i) + "), 'Color',col_array(" + num2str(i) + "))")
end
xlim([low_lim, up_lim])
%title('Koren Flux Limiter','Interpreter','latex')
ax1.TickLabelInterpreter = "latex";
xlabel('$r$','interpreter','latex')
xticks(low_lim:0.25:up_lim)
ylabel('$\Phi_{KN}$','interpreter','latex')
grid on
set(gca,'FontSize', 15)

annotation(fig1,'textarrow',[0.334 0.281],[0.309 0.262], 'Interpreter','latex',...
    'String',{'$\Phi_{KN} = 2r$'},'Color',col_array(2), 'TextLineWidth',10, 'FontSize',15);

annotation(fig1,'textarrow',[0.411 0.491],[0.666 0.626], 'Interpreter','latex',...
    'String',{'$\Phi_{KN} = \frac{2}{3} r + \frac{1}{3}$'}, 'Color',col_array(3), 'TextLineWidth',10, 'FontSize',15);

%% KOREN FUNCTION EXPORT
check_folder("pdf")
exportgraphics(fig1, exportpath() + "pdf\koren_function.pdf")

%% CFL SIMPLE VS ADVANCED
addpath('Functions')
clearvars, clc, close all
P = ParametersFunction("LE_ROY");
time_instants = logspace(0,5, 100);

[out_simple] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",0.5, "cfl_type","simple");
fprintf("Simple: %f\n", out_simple.wct_ppt)
time_func_simple = @() Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",0.5, "cfl_type","simple");
[out_advanced] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",1, "cfl_type","advanced");
time_func_advanced = @() Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","explicit", "cfl",1, "cfl_type","advanced");
fprintf("Advanced: %f\n", out_advanced.wct_ppt)

% fprintf("timeit simple: \n") % 10.0349 s
% timeit(time_func_simple)
% fprintf("timeit advanced: \n") % 6.7441 s
% timeit(time_func_advanced) % +48.7956 %

loglog(out_simple.tout, out_simple.J_Sato)
hold on
loglog(out_advanced.tout, out_advanced.J_Sato)

%% SOURCE ON VS OFF
addpath('Functions')
clearvars, clc, close all
P = ParametersFunction("LE_ROY");

time_instants = [0, logspace(0,5, 99)];

[out_source_on] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit", "source","On");
% fprintf("Source On: %f\n", out_source_on.wct + out_source_on.ppt)
[out_source_off] = Run(P, time_instants, "flux_scheme","upwind", "time_integration_scheme","implicit", "source","Off");
% fprintf("Source Off: %f\n", out_source_off.wct + out_source_off.ppt)

loglog(out_source_on.tout, out_source_on.J_Sato)
hold on
loglog(out_source_off.tout, out_source_off.J_Sato)

%%
function [] = polarization_current_plot()
    grid on
    xlabel('time ($\mathrm{s}$)', 'Interpreter','latex')
    ylabel('current density ($\mathrm{A \cdot m^{-2}}$)', 'Interpreter','latex')
    legend('Interpreter','latex')
    xticks(10.^[0 1 2 3 4 5])
    yticks(10.^[-10 -9 -8 -7 -6])
    set(gca, 'FontSize', 15)
end

function [] = computational_time_plot()
    grid on
    xlabel('simulation time ($\mathrm{s}$)', 'Interpreter','latex')
    ylabel('wall clock time ($\mathrm{s}$)', 'Interpreter','latex')
    legend('Interpreter','latex','Location','northwest')
    xticks(10.^[0 1 2 3 4 5])
    set(gca, 'FontSize', 15)
end

function [] = number_density_plot()
    grid on
    xlabel('thickness ($\mathrm{m}$)', 'Interpreter','latex')
    ylabel('number density ($\mathrm{m^{-3}}$)', 'Interpreter','latex')
    legend('Interpreter','latex')
    xlim([0 6e-5])
    set(gca, 'FontSize', 15)
end

function [] = charge_density_plot()
    grid on
    xlabel('thickness ($\mathrm{m}$)', 'Interpreter','latex')
    ylabel('charge density ($\mathrm{C \cdot m^{-3}}$)', 'Interpreter','latex')
    legend('Interpreter','latex','Location','southeast')
    set(gca, 'FontSize', 15)
end

function [plot_options_upwind, plot_options_koren] = get_plot_options_upwind_koren()
    plot_options_upwind.LineWidth = 2;
    plot_options_upwind.Color = [0, 0, 1];
    plot_options_upwind.DisplayName = 'FOU';
    
    plot_options_koren.LineWidth = 2;
    plot_options_koren.Color = [1, 0, 1];
    plot_options_koren.DisplayName = 'SOU/KL';
    plot_options_koren.LineStyle = '--';
end

function [plot_options_upwind, plot_options_koren] = extra_plot_options_upwind_koren(plot_options_upwind, plot_options_koren)
    plot_options_koren.Marker = ".";
    plot_options_koren.MarkerSize = 20;
    plot_options_upwind.Marker = ".";
    plot_options_upwind.MarkerSize = 20;
end

function [str] = my_num2str(n)
    if n == 0
        str = "";
    elseif n > 0
        str = "+" + num2str(n);
    elseif n < 0
        str = num2str(n);
    end
end
