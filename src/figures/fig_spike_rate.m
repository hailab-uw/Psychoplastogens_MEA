function fig_spike_rate(study, varargin)
%FIG_SPIKE_RATE Export spike-rate summary statistics.
%
%   FIG_SPIKE_RATE(study) loads the cached spike rates for the requested
%   study ('doi' or 'ket'), pools all baseline/treatment pairs, applies the
%   canonical channel filters, and writes CSV/JSON numeric sidecars under
%   output/fig{2,4}/stats/.
%
%   FIG_SPIKE_RATE(study, 'minRateThreshold', 1.0,                ...
%                         'ignoreSilentChannels', true,            ...
%                         'excludeMeanMultiplierOutliers', true,   ...
%                         'meanRateOutlierMultiplier', 15,         ...
%                         'channels', 1:64)
%   overrides individual options. Defaults match the proof-of-concept
%   src/figure_scripts/figures_spike_rate.m exactly.
%
% INPUTS:
%   study  -  'doi' or 'ket'.
%
% OUTPUTS:
%   CSV + JSON sidecars written to output/fig{2,4}/stats/.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels',                       cfg.channels.default);
    addParameter(p, 'ignoreSilentChannels',           true);
    addParameter(p, 'silentMode',                     cfg.silent.mode, ...
        @(s) any(strcmpi(s, {'baseline_min','both_zero','either_zero'})));
    addParameter(p, 'minRateThreshold',               cfg.silent.min_rate_spike);
    addParameter(p, 'excludeMeanMultiplierOutliers',  false);
    addParameter(p, 'meanRateOutlierMultiplier',      cfg.outlier.mean_multiplier);
    addParameter(p, 'outlierMode',                    cfg.outlier.mode, ...
        @(s) any(strcmpi(s, {'none','tukey','percentile','mean_multiplier'})));
    addParameter(p, 'outlierUpperPercentile',         cfg.outlier.upper_percentile);
    addParameter(p, 'outlierIqrFactor',               cfg.outlier.iqr_factor);
    addParameter(p, 'yScale',                         'linear', ...
        @(s) any(strcmpi(s, {'linear','log'})));
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, ~] = get_pairs_and_labels(cfg, study);
    [ratesBaseline, ratesTreatment] = load_pair_metric(pairs, opt.channels, 'spikeRates', cfg);

    if isempty(ratesBaseline)
        error('fig_spike_rate:NoData', ...
            'No valid channels found across pairs after NaN drop.');
    end

    % --- Silent-channel filter -------------------------------------------
    % Primary literature default is 'baseline_min' at 0.1 spike/s (=
    % 6 spike/min): Mossink 2021 p.2187, Brofiga 2023 p.5. Kept the
    % older 'both_zero' / 'either_zero' rules as options for sensitivity
    % supplements and for backwards compat with existing caller scripts.
    nBefore = numel(ratesBaseline);
    if opt.ignoreSilentChannels
        switch lower(opt.silentMode)
            case 'baseline_min'
                nonSilent = ratesBaseline >= opt.minRateThreshold;
            case 'both_zero'
                nonSilent = ~(ratesBaseline < opt.minRateThreshold ...
                            & ratesTreatment < opt.minRateThreshold);
            case 'either_zero'
                nonSilent = ratesBaseline  >= opt.minRateThreshold ...
                          & ratesTreatment >= opt.minRateThreshold;
            otherwise
                nonSilent = true(size(ratesBaseline));
        end
        ratesBaseline  = ratesBaseline(nonSilent);
        ratesTreatment = ratesTreatment(nonSilent);
    end
    nSilentDropped = nBefore - numel(ratesBaseline);
    fprintf('Silent filter (%s, min=%.2f/min): dropped %d of %d.\n', ...
        opt.silentMode, opt.minRateThreshold, nSilentDropped, nBefore);
    % --- Outlier filtering ------------------------------------------------
    % The primary analysis does not drop high-rate outliers. The legacy
    % mean-multiplier path is retained for sensitivity checks.
    effectiveMode = opt.outlierMode;
    if opt.excludeMeanMultiplierOutliers && strcmpi(opt.outlierMode, 'tukey')
        % Legacy behaviour: legacy flag takes precedence over default mode
        % but not over an explicit user choice.
        effectiveMode = 'mean_multiplier';
    end
    [keepMask, outlierInfo] = robust_outlier_filter(ratesBaseline, ratesTreatment, ...
        'mode',            effectiveMode, ...
        'upperPercentile', opt.outlierUpperPercentile, ...
        'multiplier',      opt.meanRateOutlierMultiplier, ...
        'iqrFactor',       opt.outlierIqrFactor);
    ratesBaseline  = ratesBaseline(keepMask);
    ratesTreatment = ratesTreatment(keepMask);
    fprintf('Outlier filter (%s): dropped %d of %d channel observation(s) [cut=%.2f].\n', ...
        outlierInfo.mode, outlierInfo.nDropped, outlierInfo.nTotal, ...
        outlierInfo.upperCutBaseline);

    if isempty(ratesBaseline)
        error('fig_spike_rate:NoData', ...
            'All channels filtered out before plotting.');
    end

    nCh         = numel(ratesBaseline);
    pctIncrease = 100 * sum(ratesTreatment > ratesBaseline) / nCh;
    pctDecrease = 100 * sum(ratesTreatment < ratesBaseline) / nCh;

    statsDir = output_path(cfg, study, 'rates', 'stats');
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end

    fprintf('fig_spike_rate(%s): n=%d, +%.1f%%, -%.1f%%\n', ...
        study, nCh, pctIncrease, pctDecrease);

    % --- Numeric sidecar for paper writing ------------------------------
    psStats = paired_stats(ratesBaseline, ratesTreatment);
    stats = struct( ...
        'study',            study, ...
        'metric',           'spike_rate', ...
        'unit',             'spikes/min', ...
        'n_channels',       nCh, ...
        'pct_increased',    pctIncrease, ...
        'pct_decreased',    pctDecrease, ...
        'median_baseline',  psStats.medianBaseline, ...
        'median_treatment', psStats.medianTreatment, ...
        'median_delta',     psStats.medianDelta, ...
        'median_pct_change',psStats.medianPctChange, ...
        'ci_pct_change_supportive', psStats.bootstrap.ciPctChange, ...
        'p_bootstrap_supportive',   psStats.bootstrap.pPctChange, ...
        'p_wilcoxon_paired_electrodes_descriptive', psStats.wilcoxon.p, ...
        'hedges_g_av',      psStats.hedgesGav);
    export_figure_stats(stats, fullfile(statsDir, ...
        sprintf('%s_spike_rate_stats', study)));
end
