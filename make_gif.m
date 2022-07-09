function make_gif(x_vect, t_vect, y_mat, filename, plottitle, plotylabel, theta_ticks)
% Make GIF from matrix

figure
gif(append('GIFs/', filename, '.gif'), 'overwrite', true)

theta_labels = strrep(string(sym(theta_ticks)), 'pi', '\pi');

% Define plot
h = plot(x_vect, y_mat(1,:));
title(plottitle)
xlabel('Angular position [rad]')
ylabel(plotylabel)
grid on;
set(gca,'XTick',theta_ticks)
set(gca,'XTickLabel', theta_labels)
ylim([min(y_mat(:)) * 1.1 max(y_mat(:)) * 1.1])     % Make sure axes are constant
set(gcf,'color','w');

% Loop through time
for i = 2:length(t_vect)
    set(h, 'Ydata', y_mat(i,:))
    gif
end

