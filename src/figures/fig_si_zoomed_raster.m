function fig_si_zoomed_raster(study, varargin)
%FIG_SI_ZOOMED_RASTER  Zoomed 20-second spike raster with burst shading.
%
%   FIG_SI_ZOOMED_RASTER(study) draws a multi-row figure for every pair
%   in STUDY ('doi' or 'ket').  Each row is one pair; the two columns are
%   baseline (left) and treatment (right).  Inside every subplot:
%
%     - Black vertical ticks mark individual spike times.
%     - Translucent coloured rectangles shade detected burst intervals
%       (burstStartTimes to burstEndTimes) behind the spike ticks, so
%       reviewers can visually resolve individual burst events.
%
%   The time window defaults to [30, 50] seconds (20 s wide); both the
%   offset and duration are configurable.
%
%   FIG_SI_ZOOMED_RASTER(study, 'Name', value, ...) options:
%
%     'pairIndex'   scalar index — show only one pair instead of all
%     'zoomSec'     width of the time window in seconds (default 20)
%     'zoomStart'   left edge of the window in seconds (default 30)
%     'tickWidth'   line width for spike ticks (default 0.35)
%     'outDir'      output directory (default output/SI/)
%     'outName'     override output filename stem
%
% INPUTS:
%   study  -  'doi' | 'ket'
%
% OUTPUTS:
%   PDF + PNG written to outDir.
%
% See also: FIG_RATE_PANEL_RASTER, LOAD_PAIR_SPIKES, LOAD_CACHE.

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study',     @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'pairIndex', [],     @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'zoomSec',   20,     @(x) isscalar(x) && x > 0);
    addParameter(p, 'zoomStart', 30,     @(x) isscalar(x) && x >= 0);
    addParameter(p, 'tickWidth', 0.35,   @(x) isscalar(x) && x > 0);
    addParameter(p, 'outDir',    output_path(cfg, '', 'si', ''));
    addParameter(p, 'outName',   '');
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    if isempty(pairs)
        error('fig_si_zoomed_raster:NoPairs', 'No pairs for study %s.', study);
    end

    % Select pairs to plot.
    if isempty(opt.pairIndex)
        pairIndices = 1:numel(pairs);
    else
        pairIndices = min(opt.pairIndex, numel(pairs));
    end
    nRows = numel(pairIndices);

    tStart = opt.zoomStart;
    tEnd   = tStart + opt.zoomSec;

    fprintf('fig_si_zoomed_raster(%s): %d pair(s), window [%.0f, %.0f] s\n', ...
        study, nRows, tStart, tEnd);

    % --- Figure layout: nRows x 2 ----------------------------------------
    figH = max(3.0, 2.4 * nRows);
    fig  = create_panel_figure(18.3, figH);
    tl   = tiledlayout(fig, nRows, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    colors = paired_plot_colors(study);
    burstAlpha = 0.30;

    for ri = 1:nRows
        pi = pairIndices(ri);
        pair = pairs(pi);

        % Load caches for burst intervals and spike times.
        bCache = load_cache(pair.baseline,  cfg);
        tCache = load_cache(pair.treatment, cfg);

        [bSpikes, tSpikes, ~, ~] = load_pair_spikes(pair, cfg.channels.default, cfg);

        % Baseline subplot.
        axB = nexttile(tl, (ri - 1) * 2 + 1);
        plot_zoomed_raster(axB, bSpikes, ...
            bCache.burstStartTimes, bCache.burstEndTimes, bCache.channelsUsed, ...
            cfg.channels.default, tStart, tEnd, opt.tickWidth, ...
            colors.baseline, burstAlpha);
        if ri == 1
            title(axB, labels.baseline, 'FontWeight', 'bold');
        end
        ylabel(axB, sprintf('Pair %d', pi), 'FontWeight', 'bold');

        % Treatment subplot.
        axT = nexttile(tl, (ri - 1) * 2 + 2);
        plot_zoomed_raster(axT, tSpikes, ...
            tCache.burstStartTimes, tCache.burstEndTimes, tCache.channelsUsed, ...
            cfg.channels.default, tStart, tEnd, opt.tickWidth, ...
            colors.treatment, burstAlpha);
        if ri == 1
            title(axT, labels.treatment, 'FontWeight', 'bold');
        end
    end

    % --- Style and save ---------------------------------------------------
    apply_nature_style(fig);

    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    if isempty(opt.outName)
        outBase = sprintf('%s_si_zoomed_raster', study);
    else
        outBase = opt.outName;
    end
    outFile = fullfile(opt.outDir, [outBase '.png']);
    save_figure(fig, outFile);
    close(fig);

    fprintf('fig_si_zoomed_raster(%s): saved %s\n', study, outBase);
end

% =========================================================================
function plot_zoomed_raster(ax, spikeTimes, burstStarts, burstEnds, ...
        chUsed, channels, tStart, tEnd, tickW, burstColor, burstAlpha)
%PLOT_ZOOMED_RASTER  Draw spike ticks and burst shading in a time window.
%
%   Burst intervals are drawn first as translucent patches; spike ticks are
%   drawn on top using the NaN-separated vertex array technique for speed.

    nCh = numel(channels);
    hold(ax, 'on');

    % --- Map cache channel indices to requested channel ordering ----------
    [~, cacheIdx] = ismember(channels, chUsed(:)');

    % --- Burst shading (patches) ------------------------------------------
    for c = 1:nCh
        ci = cacheIdx(c);
        if ci == 0; continue; end
        bst = burstStarts{ci};
        ben = burstEnds{ci};
        if isempty(bst); continue; end

        % Keep only bursts overlapping [tStart, tEnd].
        keep = bst < tEnd & ben > tStart;
        bst  = max(bst(keep), tStart);
        ben  = min(ben(keep), tEnd);
        nB   = numel(bst);
        if nB == 0; continue; end

        % Build patch vertices for all bursts on this channel at once.
        %   Each burst -> 4 corners: (bst, c-0.4), (ben, c-0.4),
        %                            (ben, c+0.4), (bst, c+0.4).
        xp = [bst(:), ben(:), ben(:), bst(:)]';  % 4 x nB
        yp = repmat([c - 0.4; c - 0.4; c + 0.4; c + 0.4], 1, nB);
        patch(ax, xp, yp, burstColor, ...
            'FaceAlpha', burstAlpha, 'EdgeColor', 'none');
    end

    % --- Spike ticks (NaN-separated line) ---------------------------------
    totalSpikes = 0;
    for c = 1:nCh
        ts = spikeTimes{c};
        totalSpikes = totalSpikes + sum(ts >= tStart & ts < tEnd);
    end

    Xall = nan(3 * totalSpikes, 1);
    Yall = nan(3 * totalSpikes, 1);
    idx  = 1;
    for c = 1:nCh
        ts = spikeTimes{c};
        ts = ts(ts >= tStart & ts < tEnd);
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
    line(ax, Xall, Yall, 'Color', 'k', 'LineWidth', tickW);

    xlim(ax, [tStart tEnd]);
    ylim(ax, [0.5, nCh + 0.5]);
    set(ax, 'YDir', 'normal');
    xlabel(ax, 'Time (s)');
    box(ax, 'on');
    hold(ax, 'off');
end
