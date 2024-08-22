function plots(U, x, t, str, eps, dx) 
    n = length(t);
    m = length(x);
    t_array = [n, ceil(3*n/4), ceil(n/2), ceil(n/4), 1];

    [x_mesh, t_mesh] = meshgrid(x,t);
    
    figure;
    mesh(x_mesh, t_mesh, U,'FaceColor','flat', 'FaceAlpha',0.6);
    title(strcat('Epsilon = ', num2str(eps), ' - Surface plot of approximations -', str));
    xlabel('S');
    ylabel('t');
    zlabel('U(S,t)');
    hold off;

    figure;
    for i = t_array
        plot(x, U(i, :), 'LineWidth',1);
        title(strcat('Epsilon = ', num2str(eps), ' - Numerical solutions at different time levels -', str));
        xlabel('x');
        ylabel('U(x,t)');
        hold on;
    end
    legend('t = 0', 't = 0.25', 't = 0.5', 't = 0.75', 't = 1');
    hold off;

    figure;
    for i = t_array
        del = zeros(1, m-1);
        for j = 1:m-1
            del(j) = (U(i,j+1) - U(i,j))/dx;
        end
        plot(x(1:m-1), del, 'LineWidth',1);
        title(strcat('Epsilon = ', num2str(eps), ' - Delta -', str));
        xlabel('x');
        ylabel('Delta');
        hold on;
    end
    legend('t = 0', 't = 0.25', 't = 0.5', 't = 0.75', 't = 1');
    hold off;

    figure;
    for i = t_array
        gamma = zeros(1, m-2);
        for j = 2:m-1
            gamma(j-1) = (U(i,j+1) - 2*U(i,j) + U(i,j-1))/(dx^2);
        end
        plot(x(2:m-1), gamma, 'LineWidth',1);
        title(strcat('Epsilon = ', num2str(eps), ' - Gamma -', str));
        xlabel('x');
        ylabel('Gamma');
        hold on;
    end
    legend('t = 0', 't = 0.25', 't = 0.5', 't = 0.75', 't = 1');
    hold off;
end