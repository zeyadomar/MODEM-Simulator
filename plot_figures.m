function plot_figures(fig_num,x_axis,y_axis,title_str,xlab,ylab,sub_ind)


figure(fig_num)
subplot(sub_ind)
plot(x_axis,y_axis)
title(title_str)
xlabel(xlab)
ylabel(ylab)





end