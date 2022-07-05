function defl = deform_plot(machine, m, A_m)
% Plot the quasi-static deflection as polar plot

% Define circle
theta = linspace(0, 2*pi, 1000);

% Increase amplitude
amplification = 10 * 10^5;
A = A_m.' * amplification;

Rs = machine.R_s;
Rs_t = Rs * ones(length(theta));

% Define modes
m0 = Rs + A(1) * cos(m(1) * theta);
m02 = Rs - A(1) * cos(m(1) * theta);
m16 = Rs + A(2) .* cos(m(2) .* theta);
mtot = Rs + sum(A(2:end) .* cos(m(2:end).' .* theta));

% Make polarplot
defl = figure;
pp.Rs = polarplot(Rs_t, 'k');
hold on
pp.m0 = polarplot(m0, 'k--');
hold on
pp.m02 = polarplot(m02, 'k--');
hold on
pp.m16 = polarplot(m16, 'r');
hold on
pp.mtot = polarplot(mtot, 'b');
grid off
Ax = gca;
Ax.RTickLabel = [];
Ax.ThetaTickLabel = [];
Ax.GridAlpha = 0.0;
hold off
legend([pp.Rs(1); pp.m0(1); pp.m16(1); pp.mtot(1)],...
    {'Stator radius','m=0', 'm='+string(m(2)), 'Total deflection'}, ...
    'Position', [0.42 0.42 0.2 0.2])
%title('Quasi-static deflection', 'FontSize', 16)

end