function plots(v, S, t, x, tau, str)
    figure;
    plot(S, v(end, :), S, v(1, :), 'LineWidth',1);
    title(strcat('Numerical solutions at final time level -', str));
    legend('t = 0', 't = T');
    xlabel('S');
    ylabel('V(S,0) and V(S,T)');
    hold off;
    saveas(gcf, strcat(str, '_V_S.png'));
    
    [S_mesh, t_mesh] = meshgrid(S,t);
    
    figure;
    mesh(S_mesh, t_mesh, v,'FaceColor','flat', 'FaceAlpha',0.6);
    title(strcat('Surface plot of approximations -', str));
    xlabel('S');
    ylabel('t');
    zlabel('V(S,t)');
    hold off;
    saveas(gcf, strcat(str, '_V_S_t.png'));

    figure;
    plot(x, v(end, :), x, v(1, :), 'LineWidth',1);
    title(strcat('Numerical solutions at final time level -', str));
    legend('tau = sigma^2/2', 'tau = 0');
    xlabel('x');
    ylabel('v(x,sigma^2/2) and v(x,0)');
    hold off;
    saveas(gcf, strcat(str, '_v_x.png'));
    
    [x_mesh, tau_mesh] = meshgrid(x,tau);
    
    figure;
    mesh(x_mesh, tau_mesh, v,'FaceColor','flat', 'FaceAlpha',0.6);
    title(strcat('Surface plot of approximations -', str));
    xlabel('x');
    ylabel('tau');
    zlabel('v(x,tau)');
    hold off;
    saveas(gcf, strcat(str, '_v_x_tau.png'));
end