function fig_si_burst_overlay_raster(study, varargin)
%FIG_SI_BURST_OVERLAY_RASTER  Full-recording raster with burst-onset markers.
%
%   FIG_SI_BURST_OVERLAY_RASTER(study) draws a two-column raster (baseline
%   on left, treatment on right) for one pair from STUDY ('doi' or 'ket').
%   Unlike the standard raster, this figure distinguishes "burst spikes"
%   from "tonic spikes":
%
%     - All spikes are drawn as thin vertical ticks in light gray.
%     - Burst onset times (burstStartTimes) are drawn as coloured vertical
%       ticks on top: red for baseline, treatment colour for treatment.
%
%   This makes it easy to see where burst events concentrate relative to
%   the full spike train.
%
%   FIG_SI_BURST_OVERLAY_RASTER(study, 'Name', value, ...) options:
%
%     'pairIndex'    which pair to show (default: middle pair)
%     'channels'     channel list (default cfg.channels.recording)
%     'tickWidth'    line width for tonic spikes (default 0.30)
%     'burstTickW'   line width for burst onset ticks (default 0.60)
%     'outDir'       output directory (default output/SI/)
%     'outName'      override output filename stem
%
% INPUTS:
%   study  -  'doi' | 'ket'
%
% OUTPUTS:
%   PDF + PNG written to outDir.
%
% See also: FIG_RATE_PANEL_RASTER, FIG_SI_ZOOMED_RASTER, LOAD_CACHE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study',      @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'pairIndex',  [],    @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'channels',   cfg.channels.recording);
    addParameter(p, 'tickWidth',  0.30,  @(x) isscalar(x) && x > 0);
    addParameter(p, 'burstTickW', 0.60,  @(x) isscalar(x) && x > 0);
    addParameter(p, 'outDir',     output_path(cfg, '', 'si', ''));
    addParameter(p, 'outName',    '');
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    if isempty(pairs)
        error('fig_si_burst_overlay_raster:NoPairs', 'No pairs for study %s.', study);
    end
    if isempty(opt.pairIndex)
        pairIdx = max(1, round(numel(pairs) / 2));
    else
        pairIdx = min(opt.pairIndex, numel(pairs));
    end
    pair = pairs(pairIdx);

    fprintf(['fig_si_burst_overlay_raster(%s): pair %d of %d\n' ...
             '  baseline:  %s\n  treatment: %s\n'], ...
        study, pairIdx, numel(pairs), pair.baseline, pair.treatment);

    % --- Load spike times and burst onset times ---------------------------
    [bSpikes, tSpikes, bMeta, tMeta] = load_pair_spikes(pair, opt.channels, cfg);
    bCache = load_cache(pair.baseline,  cfg);
    tCache = load_cache(pair.treatment, cfg);

    bBurstOnsets = align_burst_onsets(bCache.burstStartTimes, bCache.channelsUsed, opt.channels);
    tBurstOnsets = align_burst_onsets(tCache.burstStartTimes, tCache.channelsUsed, opt.channels);

    colors = paired_plot_colors(study);

    % --- Figure: 1 x 2 layout (baseline | treatment) ---------------------
    fig = create_panel_figure(18.3, 5.0);
    tl  = tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('Burst overlay raster  (pair %d)', pairIdx), ...
        'FontWeight', 'bold', 'Interpreter', 'none');

    burstColorB = [0.80 0.20 0.20];   % red for baseline burst onsets
    burstColorT = colors.treatment;    % study treatment colour

    axB = nexttile(tl, 1);
    plot_overlay_raster(axB, bSpikes, bBurstOnsets, bMeta.durationSec, ...
        labels.baseline, opt.tickWidth, opt.burstTickW, burstColorB);

    axT = nexttile(tl, 2);
    plot_overlay_raster(axT, tSpikes, tBurstOnsets, tMeta.durationSec, ...
        labels.treatment, opt.tickWidth, opt.burstTickW, burstColorT);

    % --- Style and save ---------------------------------------------------
    apply_nature_style(fig);

    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    if isempty(opt.outName)
        outBase = sprintf('%s_si_burst_overlay', study);
    else
        outBase = opt.outName;
    end
    outFile = fullfile(opt.outDir, [outBase '.png']);
    save_figure(fig, outFile);
    close(fig);

    fprintf('fig_si_burst_overlay_raster(%s): saved %s\n', study, outBase);
end

% =========================================================================
function onsets = align_burst_onsets(burstStartTimes, chUsed, channels)
%ALIGN_BURST_ONSETS  Return burst onset cell array aligned to CHANNELS.
    nCh    = numel(channels);
    onsets = repmat({zeros(0, 1)}, 1, nCh);
    [~, idx] = ismember(channels, chUsed(:)');
    for k = 1:nCh
        if idx(k) > 0
            v = burstStartTimes{idx(k)};
            onsets{k} = v(:);
        end
    end
end

% =========================================================================
function plot_overlay_raster(ax, spikeTimes, burstOnsets, durSec, ...
        ttl, tickW, burstTickW, burstColor)
%PLOT_OVERLAY_RASTER  Tonic spikes in gray, burst onsets in colour.

    nCh = numel(spikeTimes);
    hold(ax, 'on');

    % --- Tonic spikes: all spikes in light gray ---------------------------
    totalSpikes = sum(cellfun(@numel, spikeTimes));
    Xall = nan(3 * totalSpikes, 1);
    Yall = nan(3 * totalSpikes, 1);
    idx  = 1;
    for c = 1:nCh
        ts = spikeTimes{c};
        n  = numel(ts);
        if n == 0; continue; end
        k = idx:(idx + 3*n - 1);
        Xall(k(1:3:end)) = ts(:);
        Xall(k(2:3:end)) = ts(:);
        Xall(k(3:3:end)) = NaN;
        Yall(k(1:3:end)) = c - 0.4;
        Yall(k(2:3:end)) = c + 0.4;
        Yall(k(3:3:end)) = NaN;
        idx = idx + 3*n;
    end
    line(ax, Xall, Yall, 'Color', [0.78 0.78 0.78], 'LineWidth', tickW);

    % --- Burst onsets: coloured ticks on top ------------------------------
    totalOnsets = sum(cellfun(@numel, burstOnsets));
    Xb = nan(3 * totalOnsets, 1);
    Yb = nan(3 * totalOnsets, 1);
    idx = 1;
    for c = 1:nCh
        bo = burstOnsets{c};
        n  = numel(bo);
        if n == 0; continue; end
        k = idx:(idx + 3*n - 1);
        Xb(k(1:3:end)) = bo(:);
        Xb(k(2:3:end)) = bo(:);
        Xb(k(3:3:end)) = NaN;
        Yb(k(1:3:end)) = c - 0.4;
        Yb(k(2:3:end)) = c + 0.4;
        Yb(k(3:3:end)) = NaN;
        idx = idx + 3*n;
    end
    line(ax, Xb, Yb, 'Color', burstColor, 'LineWidth', burstTickW);

    xlim(ax, [0 max(durSec, 1)]);
    ylim(ax, [0.5, nCh + 0.5]);
    set(ax, 'YDir', 'normal');
    xlabel(ax, 'Time (s)', 'FontWeight', 'bold');
    ylabel(ax, 'Channel',  'FontWeight', 'bold');
    title(ax, ttl, 'FontWeight', 'bold');
    box(ax, 'on');
    hold(ax, 'off');
end
