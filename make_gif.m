function make_gif(x_vect, t_vect, y_mat, filename)
%Make GIF from matrix

figure
gif(append('GIFs/', filename, '.gif'), 'overwrite', true)

h = plot(x_vect, y_mat(1,:));
ylim([min(y_mat(:)) * 1.1 max(y_mat(:)) * 1.1])
for i = 2:length(t_vect)
    set(h, 'Ydata', y_mat(i,:))
    gif
end

