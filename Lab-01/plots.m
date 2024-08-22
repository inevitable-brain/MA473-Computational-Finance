function plots(x,t,approximations,actual_f,str)
    figure;
    plot(x,approximations(length(t),:));
    hold on;
    plot(x,actual_f(x,1));
    hold off;
    grid on;
    legend('Numerical','Exact');
    xlabel('X');
    ylabel('u(x,t)');
    title(str);


    [x1, t1] = meshgrid(x, t);
    u = actual_f(x1, t1);
    
    figure;
    
    % Plot the approximations surface
    surf(x, t, approximations, 'FaceColor', 'y', 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    hold on;
    % Plot the actual function surface
    surf(x1, t1, u, 'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    
    xlabel('X');
    ylabel('t');
    zlabel('Exact and Numerical Solutions');
    
    legend('Approximations', 'Actual');
    title('Exact and Numerical plots');
    
    hold off;

end