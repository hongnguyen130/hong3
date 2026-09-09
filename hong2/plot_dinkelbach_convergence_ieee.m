function stats = plot_dinkelbach_convergence_ieee(residual_hist, tol, savePrefix)
% Plot Dinkelbach convergence in an IEEE-style manner
%
% INPUTS:
%   residual_hist : Ns x maxIter matrix
%                   each row = one Monte Carlo run
%                   values should be abs(SSR - r*Ptot)
%                   use NaN for iterations not reached
%
%   tol           : convergence threshold, e.g. 1e-4 or 1e-5
%   savePrefix    : string for saving figures, e.g. 'fig/dinkelbach'
%                   if empty, figures are not saved
%
% OUTPUT:
%   stats         : structure containing summary statistics
%
% Example:
%   stats = plot_dinkelbach_convergence_ieee(residual_hist, 1e-4, 'dink_conv');

    if nargin < 2 || isempty(tol)
        tol = 1e-4;
    end
    if nargin < 3
        savePrefix = '';
    end

    % -----------------------------
    % Basic sizes
    % -----------------------------
    [Ns, maxIter] = size(residual_hist);

    % -----------------------------
    % Clean data
    % -----------------------------
    R = residual_hist;

    % Replace non-positive values by eps for log plotting
    R(R <= 0) = eps;

    % -----------------------------
    % Detect convergence iteration
    % k_conv(i) = first iteration where residual <= tol
    % if never converged, set to maxIter
    % -----------------------------
    k_conv = maxIter * ones(Ns,1);
    is_conv = false(Ns,1);

    for i = 1:Ns
        idx = find(~isnan(R(i,:)) & R(i,:) <= tol, 1, 'first');
        if ~isempty(idx)
            k_conv(i) = idx;
            is_conv(i) = true;
        else
            % if run has trailing NaNs, use last valid iteration
            lastValid = find(~isnan(R(i,:)), 1, 'last');
            if ~isempty(lastValid)
                k_conv(i) = lastValid;
            end
        end
    end

    % -----------------------------
    % Padding after convergence
    % This is the key step to avoid fake zig-zag in median curves
    % -----------------------------
    R_pad = R;
    for i = 1:Ns
        lastValid = find(~isnan(R_pad(i,:)), 1, 'last');
        if isempty(lastValid)
            continue;
        end

        kc = k_conv(i);
        kc = min(kc, lastValid);

        % Fill NaNs before lastValid if needed (rare case)
        for k = 2:lastValid
            if isnan(R_pad(i,k))
                R_pad(i,k) = R_pad(i,k-1);
            end
        end

        % Pad after convergence using converged residual
        if kc < maxIter
            R_pad(i, kc+1:maxIter) = R_pad(i, kc);
        end

        % If there are NaNs after lastValid, also pad them
        if lastValid < maxIter
            R_pad(i, lastValid+1:maxIter) = R_pad(i, lastValid);
        end
    end

    % Final safety
    R_pad(isnan(R_pad)) = eps;
    R_pad(R_pad <= 0) = eps;

    % -----------------------------
    % Statistics across runs
    % -----------------------------
    medR = median(R_pad, 1);
    p10  = prctile(R_pad, 10, 1);
    p90  = prctile(R_pad, 90, 1);

    x = 1:maxIter;

    % -----------------------------
    % Figure 1: Median + percentile band
    % -----------------------------
    fig1 = figure('Color', 'w', 'Position', [100 100 760 560]);
    hold on; box on; grid on;

    % percentile band
    Xfill = [x, fliplr(x)];
    Yfill = [p10, fliplr(p90)];
    hBand = fill(Xfill, Yfill, [0.85 0.90 1.00], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.35);

    % median
    hMed = semilogy(x, medR, '-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'w');

    % percentile lines
    hP10 = semilogy(x, p10, '--', 'LineWidth', 1.2);
    hP90 = semilogy(x, p90, '--', 'LineWidth', 1.2);

    % tolerance line
    hTol = line(xlim,[tol tol],'LineStyle', ':', 'LineWidth', 1.2);

    xlabel('Dinkelbach iteration', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel('|SSR - rP_{tot}|', 'FontName', 'Times New Roman', 'FontSize', 13);
    title(sprintf('Dinkelbach convergence over %d Monte Carlo runs', Ns), ...
        'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    legend([hMed, hP10, hP90, hBand, hTol], ...
        {'Median', '10th percentile', '90th percentile', '10-90 percentile band', 'Tolerance'}, ...
        'Location', 'southwest', 'FontName', 'Times New Roman', 'FontSize', 11);

    set(gca, 'FontName', 'Times New Roman', ...
             'FontSize', 12, ...
             'LineWidth', 1.0, ...
             'GridAlpha', 0.20, ...
             'MinorGridAlpha', 0.15, ...
             'YMinorGrid', 'on', ...
             'XMinorGrid', 'off');

    axis tight;

    % -----------------------------
    % Figure 2: Median only (clean paper version)
    % -----------------------------
    fig2 = figure('Color', 'w', 'Position', [140 120 700 520]);
    semilogy(x, medR, '-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 6, ...
        'MarkerFaceColor', 'w');
    hold on; box on; grid on;
    line(xlim, [tol tol], 'LineStyle',':', 'LineWidth', 1.2);

    xlabel('Dinkelbach iteration', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel('|SSR - rP_{tot}|', 'FontName', 'Times New Roman', 'FontSize', 13);
    title('Dinkelbach convergence (median trajectory)', ...
        'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    set(gca, 'FontName', 'Times New Roman', ...
             'FontSize', 12, ...
             'LineWidth', 1.0, ...
             'GridAlpha', 0.20, ...
             'MinorGridAlpha', 0.15, ...
             'YMinorGrid', 'on');

    axis tight;

    % -----------------------------
    % Figure 3: CDF of convergence iteration
    % -----------------------------
    fig3 = figure('Color', 'w', 'Position', [180 140 700 520]);
    [f, xx] = ecdf(k_conv);
    stairs(xx, f, 'LineWidth', 1.8);
    hold on; box on; grid on;

    xlabel('Convergence iteration', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel('Empirical CDF', 'FontName', 'Times New Roman', 'FontSize', 13);
    title('CDF of convergence iteration', ...
        'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    set(gca, 'FontName', 'Times New Roman', ...
             'FontSize', 12, ...
             'LineWidth', 1.0, ...
             'GridAlpha', 0.20);

    xlim([1 maxIter]);
    ylim([0 1]);

    % -----------------------------
    % Figure 4: Histogram of convergence iteration
    % -----------------------------
    fig4 = figure('Color', 'w', 'Position', [220 160 700 520]);
    histogram(k_conv, 'BinMethod', 'integers');
    hold on; box on; grid on;

    xlabel('Convergence iteration', 'FontName', 'Times New Roman', 'FontSize', 13);
    ylabel('Number of runs', 'FontName', 'Times New Roman', 'FontSize', 13);
    title('Histogram of convergence iteration', ...
        'FontName', 'Times New Roman', 'FontSize', 13, 'FontWeight', 'bold');

    set(gca, 'FontName', 'Times New Roman', ...
             'FontSize', 12, ...
             'LineWidth', 1.0, ...
             'GridAlpha', 0.20);

    xlim([1 maxIter]);

    % -----------------------------
    % Summary statistics
    % -----------------------------
    stats = struct();
    stats.Ns = Ns;
    stats.maxIter = maxIter;
    stats.tol = tol;
    stats.k_conv = k_conv;
    stats.is_conv = is_conv;
    stats.convRate = mean(is_conv);
    stats.meanConvIter = mean(k_conv);
    stats.medianConvIter = median(k_conv);
    stats.medR = medR;
    stats.p10 = p10;
    stats.p90 = p90;
    stats.R_pad = R_pad;

    fprintf('\n===== Dinkelbach convergence summary =====\n');
    fprintf('Number of Monte Carlo runs     : %d\n', Ns);
    fprintf('Maximum iterations             : %d\n', maxIter);
    fprintf('Tolerance                      : %.2e\n', tol);
    fprintf('Convergence rate               : %.2f %%\n', 100*stats.convRate);
    fprintf('Mean convergence iteration     : %.2f\n', stats.meanConvIter);
    fprintf('Median convergence iteration   : %.2f\n', stats.medianConvIter);

    % -----------------------------
    % Save figures if requested
    % -----------------------------
    if ~isempty(savePrefix)
        exportgraphics(fig1, [savePrefix '_band.pdf'], 'ContentType', 'vector');
        exportgraphics(fig2, [savePrefix '_median.pdf'], 'ContentType', 'vector');
        exportgraphics(fig3, [savePrefix '_cdf.pdf'], 'ContentType', 'vector');
        exportgraphics(fig4, [savePrefix '_hist.pdf'], 'ContentType', 'vector');

        exportgraphics(fig1, [savePrefix '_band.png'], 'Resolution', 400);
        exportgraphics(fig2, [savePrefix '_median.png'], 'Resolution', 400);
        exportgraphics(fig3, [savePrefix '_cdf.png'], 'Resolution', 400);
        exportgraphics(fig4, [savePrefix '_hist.png'], 'Resolution', 400);
    end
end