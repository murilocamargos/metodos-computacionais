function grafico( x, u )
    n = size(u, 2);
    plot(repmat(x, length(n), 1)', u);

    legends = {};
    for i = 1:n
        legends{i} = ['n = ' num2str(i)];
    end
    legend(legends, 'Location', 'southwest');
    grid on;
end

