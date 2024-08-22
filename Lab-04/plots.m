function plots(U, x, t, x_max, str, dx) 
    idx = find(x > x_max);
    x = x(1:idx(1));
    U = U(:, 1:idx(1));
    
    n = length(t);
    m = length(x);
    t_array = [n, ceil(3*n/4), ceil(n/2), ceil(n/4), 1];

    [x_mesh, t_mesh] = meshgrid(x,t);
    
    figure;
    mesh(x_mesh, t_mesh, U,'FaceColor','flat', 'FaceAlpha',0.6);
    title(strcat('Surface plot of approximations -', str));
    xlabel('S');
    ylabel('t');
    zlabel('V(S,t)');
    hold off;
    saveas(gcf, strcat('Surf_plot_', str, '.png'));

    figure;
    for i = t_array
        plot(x, U(i, :), 'LineWidth',1);
        title(strcat('Numerical solutions at different time levels -', str));
        xlabel('S');
        ylabel('t');
        hold on;
    end
    legend('t = 0', 't = 0.25', 't = 0.5', 't = 0.75', 't = 1');
    hold off;
    saveas(gcf, strcat('Num_soln_', str, '.png'));

    figure;
    for i = t_array
        del = zeros(1, m-1);
        for j = 1:m-1
            del(j) = (U(i,j+1) - U(i,j))/dx;
        end
        plot(x(1:m-1), del, 'LineWidth',1);
        title(strcat('Delta -', str));
        xlabel('x');
        ylabel('Delta');
        hold on;
    end
    legend('t = 0', 't = 0.25', 't = 0.5', 't = 0.75', 't = 1');
    hold off;
    saveas(gcf, strcat('Delta', str, '.png'));

    figure;
    for i = t_array
        gamma = zeros(1, m-2);
        for j = 2:m-1
            gamma(j-1) = (U(i,j+1) - 2*U(i,j) + U(i,j-1))/(dx^2);
        end
        plot(x(2:m-1), gamma, 'LineWidth',1);
        title(strcat('Gamma -', str));
        xlabel('x');
        ylabel('Gamma');
        hold on;
    end
    legend('t = 0', 't = 0.25', 't = 0.5', 't = 0.75', 't = 1');
    hold off;
    saveas(gcf, strcat('Gamma', str, '.png'));
end