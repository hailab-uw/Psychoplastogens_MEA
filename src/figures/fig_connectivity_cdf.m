function fig_connectivity_cdf(study, varargin)
%FIG_CONNECTIVITY_CDF  Edge-weight distribution shift panel (Panel a).
%
%   FIG_CONNECTIVITY_CDF(study) generates a CDF of all pairwise cross-
%   correlation peak values, overlaying baseline vs treatment, pooled
%   across all recording pairs. This replaces the exemplar cross-
%   correlograms (now in SI) and shows the GLOBAL shift in edge weights
%   without cherry-picking.
%
%   FIG_CONNECTIVITY_CDF(study, 'results', R) uses pre-computed results
%   from run_connectivity, skipping the expensive recomputation.
%
% See also: RUN_CONNECTIVITY, FIG_CONNECTIVITY_SUMMARY, RUN_FIGURES.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'results',       [],  @(x) isempty(x) || isstruct(x));
    addParameter(p, 'channels',      cfg.channels.recording);
    addParameter(p, 'binMs',         cfg.connectivity.bin_ms);
    addParameter(p, 'maxLagMs',      cfg.connectivity.max_lag_ms);
    addParameter(p, 'normalization', cfg.connectivity.normalization);
    addParameter(p, 'edgeThreshold', cfg.connectivity.edge_threshold);
    addParameter(p, 'edgeDensity',   cfg.connectivity.edge_density);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [~, labels] = get_pairs_and_labels(cfg, study);
    colors = paired_plot_colors(study);

    % --- Get connectivity results -----------------------------------------
    if isempty(opt.results)
        results = run_connectivity(study, ...
            'channels',      opt.channels, ...
            'binMs',         opt.binMs, ...
            'maxLagMs',      opt.maxLagMs, ...
            'normalization', opt.normalization, ...
            'edgeThreshold', opt.edgeThreshold, ...
            'edgeDensity',   opt.edgeDensity);
    else
        results = opt.results;
    end

    nPairs = numel(results);
    nCh    = size(results(1).baseline.adjacency, 1);
    upper  = triu(true(nCh), 1);

    % --- Collect all upper-triangle weights per condition ------------------
    allBaseline  = [];
    allTreatment = [];
    for k = 1:nPairs
        bVals = results(k).baseline.adjacency(upper);
        tVals = results(k).treatment.adjacency(upper);
        allBaseline  = [allBaseline;  bVals(~isnan(bVals))]; %#ok<AGROW>
        allTreatment = [allTreatment; tVals(~isnan(tVals))]; %#ok<AGROW>
    end

    % --- Compute CDFs -----------------------------------------------------
    xRange = linspace(min([allBaseline; allTreatment]), ...
                      max([allBaseline; allTreatment]), 500);
    cdfB = arrayfun(@(x) mean(allBaseline  <= x), xRange);
    cdfT = arrayfun(@(x) mean(allTreatment <= x), xRange);

    medB = median(allBaseline);
    medT = median(allTreatment);

    % 10% density threshold (90th percentile of pooled baseline)
    thresh90 = prctile(allBaseline, 100 * (1 - opt.edgeDensity));

    % --- Plot CDF ---------------------------------------------------------
    fig = create_panel_figure(7.0, 5.5);
    ax  = axes(fig);
    hold(ax, 'on');

    plot(ax, xRange, cdfB, '-', 'Color', colors.baseline,  'LineWidth', 1.4);
    plot(ax, xRange, cdfT, '-', 'Color', colors.treatment, 'LineWidth', 1.4);

    % Median dashed lines
    yl = [0 1];
    plot(ax, [medB medB], yl, '--', 'Color', colors.baseline,  'LineWidth', 0.6);
    plot(ax, [medT medT], yl, '--', 'Color', colors.treatment, 'LineWidth', 0.6);

    % 10% density threshold marker
    plot(ax, [thresh90 thresh90], yl, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6);
    text(ax, thresh90, 0.12, sprintf(' 10%% density\n threshold'), ...
        'FontSize', 6, 'FontName', 'Arial', 'Color', [0.4 0.4 0.4]);

    % Annotation
    nEdges = sum(~isnan(results(1).baseline.adjacency(upper)));
    text(ax, 0.97, 0.28, ...
        sprintf('n_{edges} = %d per pair\nn_{pairs} = %d\n\\Delta median = %.3f', ...
            nEdges, nPairs, medT - medB), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'FontSize', 7, 'FontName', 'Arial');

    xlabel(ax, 'Peak cross-correlation (z)');
    ylabel(ax, 'Cumulative proportion');
    title(ax, 'Edge-weight distribution');
    legend(ax, {labels.baseline, labels.treatment}, ...
        'Location', 'southeast', 'Box', 'off');
    set(ax, 'TickDir', 'out', 'Box', 'off');
    ylim(ax, yl);
    hold(ax, 'off');

    apply_nature_style(fig);

    % --- Save -------------------------------------------------------------
    outDir = output_path(cfg, study, 'connectivity', '');
    if ~exist(outDir, 'dir'); mkdir(outDir); end
    save_figure(fig, fullfile(outDir, [study '_edge_weight_cdf']));
    close(fig);

    % --- Numeric sidecar --------------------------------------------------
    statsDir = output_path(cfg, study, 'connectivity', 'stats');
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end
    sidecar = struct( ...
        'study',                    study, ...
        'n_pairs',                  nPairs, ...
        'n_edges_per_pair',         nEdges, ...
        'total_edges',              numel(allBaseline), ...
        'mean_baseline',            mean(allBaseline), ...
        'mean_treatment',           mean(allTreatment), ...
        'median_baseline',          medB, ...
        'median_treatment',         medT, ...
        'delta_mean',               mean(allTreatment) - mean(allBaseline), ...
        'delta_median',             medT - medB, ...
        'density_threshold_baseline', thresh90);
    export_figure_stats(sidecar, ...
        fullfile(statsDir, [study '_edge_weight_cdf_stats']));

    fprintf('fig_connectivity_cdf(%s): saved to %s\n', study, outDir);
end
