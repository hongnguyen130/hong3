function plot_dinkelbach_onefig(residual_hist, tol)

    if nargin < 2
        tol = 1e-4;
    end

    [Ns, maxIter] = size(residual_hist);
    R = residual_hist;
    R(R <= 0) = eps;

    % find convergence iteration and pad
    for i = 1:Ns
        idx = find(~isnan(R(i,:)) & R(i,:) <= tol, 1, 'first');
        if isempty(idx)
            idx = find(~isnan(R(i,:)), 1, 'last');
        end
        if isempty(idx)
            continue;
        end

        for k = 2:maxIter
            if isnan(R(i,k))
                R(i,k) = R(i,k-1);
            end
        end

        if idx < maxIter
            R(i, idx+1:maxIter) = R(i, idx);
        end
    end

    R(isnan(R)) = eps;

    medR = median(R,1);
    p10  = prctile(R,10,1);
    p90  = prctile(R,90,1);
    x = 1:maxIter;

    figure('Color','w','Position',[100 100 760 560]);
    hold on; box on; grid on;

    fill([x fliplr(x)], [p10 fliplr(p90)], [0.86 0.90 1.00], ...
         'EdgeColor','none', 'FaceAlpha',0.35);

    semilogy(x, medR, '-o', 'LineWidth',1.8, 'MarkerSize',6, 'MarkerFaceColor','w');
    semilogy(x, p10, '--', 'LineWidth',1.2);
    semilogy(x, p90, '--', 'LineWidth',1.2);
    yline(tol, ':', 'LineWidth',1.2);

    xlabel('Dinkelbach iteration', 'FontName','Times New Roman', 'FontSize',13);
    ylabel('|SSR - rP_{tot}|', 'FontName','Times New Roman', 'FontSize',13);
    title(sprintf('Dinkelbach convergence over %d Monte Carlo runs', Ns), ...
          'FontName','Times New Roman', 'FontSize',13, 'FontWeight','bold');

    legend('10-90 percentile band', 'Median', '10th percentile', ...
           '90th percentile', 'Tolerance', ...
           'Location','southwest', ...
           'FontName','Times New Roman', 'FontSize',11);

    set(gca, 'FontName','Times New Roman', ...
             'FontSize',12, ...
             'LineWidth',1.0, ...
             'GridAlpha',0.2, ...
             'YMinorGrid','on');

    axis tight;
end