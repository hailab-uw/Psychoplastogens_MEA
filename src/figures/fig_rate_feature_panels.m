function fig_rate_feature_panels(study, varargin)
%FIG_RATE_FEATURE_PANELS MFR, MBR, IBI, and participation panels.
%
%   FIG_RATE_FEATURE_PANELS(study) writes publication panels for the
%   requested study:
%       <study>_rate_change_summary.png
%       <study>_mfr_paired_density.png
%       <study>_mbr_paired_density.png
%       <study>_ibi_paired_density.png
%       <study>_bursting_electrodes.png

    cfg = project_config();

    p = inputParser;
    addRequired(p,  'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'channels', cfg.channels.recording);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, labels] = get_pairs_and_labels(cfg, study);
    if isempty(pairs)
        error('fig_rate_feature_panels:NoPairs', 'No pairs for study %s.', study);
    end

    panelDir = output_path(cfg, study, 'rates', '');
    statsDir = output_path(cfg, study, 'rates', 'stats');
    if ~exist(panelDir, 'dir'); mkdir(panelDir); end
    if ~exist(statsDir, 'dir'); mkdir(statsDir); end

    style = figure_style_config();
    colors = rate_panel_colors(study, style);
    rateData = collect_rate_feature_data(pairs, opt.channels, cfg);

    plot_rate_change_summary(rateData, colors, style, ...
        fullfile(panelDir, sprintf('%s_rate_change_summary.png', study)));
    plot_paired_rate_density(rateData.spikeBaseline, rateData.spikeTreatment, ...
        'mfr', labels.treatment, colors, style, ...
        fullfile(panelDir, sprintf('%s_mfr_paired_density.png', study)));
    plot_paired_rate_density(rateData.burstBaseline, rateData.burstTreatment, ...
        'mbr', labels.treatment, colors, style, ...
        fullfile(panelDir, sprintf('%s_mbr_paired_density.png', study)));
    plot_paired_rate_density(rateData.ibiBaseline, rateData.ibiTreatment, ...
        'ibi', labels.treatment, colors, style, ...
        fullfile(panelDir, sprintf('%s_ibi_paired_density.png', study)));
    plot_bursting_electrodes(rateData, labels.treatment, colors, style, ...
        fullfile(panelDir, sprintf('%s_bursting_electrodes.png', study)));

    stats = summarize_rate_features(rateData, study);
    export_figure_stats(stats, fullfile(statsDir, sprintf('%s_rate_feature_stats', study)));
end

% =========================================================================
function data = collect_rate_feature_data(pairs, channels, cfg)
    nPairs = numel(pairs);
    data.mfrPct = nan(nPairs, 1);
    data.mbrPct = nan(nPairs, 1);
    data.activeElectrodesBaseline = nan(nPairs, 1);
    data.activeElectrodesTreatment = nan(nPairs, 1);
    data.burstingElectrodesBaseline = nan(nPairs, 1);
    data.burstingElectrodesTreatment = nan(nPairs, 1);
    data.spikeBaseline = [];
    data.spikeTreatment = [];
    data.burstBaseline = [];
    data.burstTreatment = [];
    data.ibiBaseline = [];
    data.ibiTreatment = [];

    for k = 1:nPairs
        [loadedB, loadedT] = load_pair_cache(pairs(k), cfg);
        [spikeB, spikeT] = align_metric(loadedB, loadedT, channels, 'spikeRates');
        [burstB, burstT] = align_metric(loadedB, loadedT, channels, 'burstRates');
        [ibiB, ibiT] = align_burst_feature(loadedB, loadedT, channels, 'ibi');

        mfrUse = spikeB >= cfg.silent.min_rate_spike;
        mbrUse = burstB > 0;

        data.mfrPct(k) = median(percent_change(spikeT(mfrUse), spikeB(mfrUse)), 'omitnan');
        data.mbrPct(k) = median(percent_change(burstT(mbrUse), burstB(mbrUse)), 'omitnan');

        data.activeElectrodesBaseline(k) = sum(spikeB >= cfg.silent.min_rate_spike);
        data.activeElectrodesTreatment(k) = sum(spikeT >= cfg.silent.min_rate_spike);
        data.burstingElectrodesBaseline(k) = sum(burstB > 0);
        data.burstingElectrodesTreatment(k) = sum(burstT > 0);

        data.spikeBaseline = [data.spikeBaseline; spikeB(:)]; %#ok<AGROW>
        data.spikeTreatment = [data.spikeTreatment; spikeT(:)]; %#ok<AGROW>
        data.burstBaseline = [data.burstBaseline; burstB(:)]; %#ok<AGROW>
        data.burstTreatment = [data.burstTreatment; burstT(:)]; %#ok<AGROW>
        data.ibiBaseline = [data.ibiBaseline; ibiB(:)]; %#ok<AGROW>
        data.ibiTreatment = [data.ibiTreatment; ibiT(:)]; %#ok<AGROW>
    end
end

function [b, t] = align_metric(loadedB, loadedT, channels, metric)
    valuesB = loadedB.(metric)(:);
    valuesT = loadedT.(metric)(:);
    chB = loadedB.channelsUsed(:)';
    chT = loadedT.channelsUsed(:)';

    [~, idxB] = ismember(channels, chB);
    [~, idxT] = ismember(channels, chT);
    valid = (idxB > 0) & (idxT > 0);
    b = valuesB(idxB(valid));
    t = valuesT(idxT(valid));

    ok = isfinite(b) & isfinite(t);
    b = b(ok);
    t = t(ok);
end

function [b, t] = align_burst_feature(loadedB, loadedT, channels, feature)
    chB = loadedB.channelsUsed(:)';
    chT = loadedT.channelsUsed(:)';

    [~, idxB] = ismember(channels, chB);
    [~, idxT] = ismember(channels, chT);
    valid = (idxB > 0) & (idxT > 0);
    idxB = idxB(valid);
    idxT = idxT(valid);

    b = nan(numel(idxB), 1);
    t = nan(numel(idxT), 1);
    for i = 1:numel(idxB)
        b(i) = channel_burst_feature(loadedB, idxB(i), feature);
        t(i) = channel_burst_feature(loadedT, idxT(i), feature);
    end
end

function value = channel_burst_feature(cacheStruct, idx, feature)
    starts = cacheStruct.burstStartTimes{idx}(:);
    switch feature
        case 'ibi'
            starts = starts(isfinite(starts));
            if numel(starts) < 2
                value = NaN;
                return;
            end
            intervalsSec = diff(sort(starts));
            intervalsSec = intervalsSec(isfinite(intervalsSec) & intervalsSec > 0);
            value = median(intervalsSec, 'omitnan');
        otherwise
            error('channel_burst_feature:UnknownFeature', ...
                'Unknown burst feature: %s', feature);
    end
end

function pct = percent_change(treatment, baseline)
    pct = nan(size(baseline));
    ok = baseline > 0;
    pct(ok) = 100 * (treatment(ok) - baseline(ok)) ./ baseline(ok);
end

function colors = rate_panel_colors(study, style)
    colors.baseline = style.baselineColor;
    colors.burst = style.burstColor;
    colors.increase = style.increaseColor;
    colors.decrease = style.decreaseColor;
    colors.silenced = style.silencedColor;
    colors.gained = style.gainedColor;
    colors.neutralLine = style.neutralLineColor;

    switch lower(study)
        case 'doi'
            colors.treatment = style.doiColor;
        case 'ket'
            colors.treatment = style.ketColor;
        otherwise
            error('rate_panel_colors:UnknownStudy', 'study must be doi or ket.');
    end
end

function plot_rate_change_summary(data, colors, style, outFile)
    fig = create_panel_figure(style.ratePanelWidthCm, style.ratePanelHeightCm);
    ax = axes(fig);
    hold(ax, 'on');

    x = (1:numel(data.mfrPct))';
    xMfr = x;
    xMbr = x;
    plot(ax, xMfr, data.mfrPct, '-', ...
        'Color', style.neutralLineColor, 'LineWidth', style.pairedLineWidth);
    plot(ax, xMbr, data.mbrPct, '-', ...
        'Color', style.neutralLineColor, 'LineWidth', style.pairedLineWidth);
    hMfr = scatter(ax, xMfr, data.mfrPct, style.summaryMarkerSize, colors.treatment, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', style.markerEdgeWidth);
    hMbr = scatter(ax, xMbr, data.mbrPct, style.summaryMarkerSize, colors.burst, 's', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', style.markerEdgeWidth);
    yline(ax, 0, '-', 'Color', style.zeroLineColor, 'LineWidth', style.referenceLineWidth, 'HandleVisibility', 'off');

    xlim(ax, [0.5, numel(x) + 0.5]);
    set(ax, 'XTick', x, 'XTickLabel', compose('%d', x));
    ylim(ax, padded_limits([data.mfrPct(:); data.mbrPct(:)], [-110 330]));
    xlabel(ax, 'Well');
    ylabel(ax, 'Median change (%)');
    legend(ax, [hMfr, hMbr], {'MFR', 'MBR'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', style.legendFontSize);
    style_axis(ax, style);
    save_and_close(fig, outFile, style);
end

function plot_paired_rate_density(baseline, treatment, metric, treatmentLabel, colors, style, outFile)
    switch metric
        case 'mfr'
            keep = baseline >= 6 | treatment >= 6;
            yLabel = 'log(MFR)';
        case 'mbr'
            keep = baseline > 0 | treatment > 0;
            yLabel = 'log(MBR)';
        case 'ibi'
            keep = isfinite(baseline) & isfinite(treatment) & baseline > 0 & treatment > 0;
            yLabel = 'log(IBI)';
        otherwise
            error('plot_paired_rate_density:UnknownMetric', ...
                'Unknown metric: %s', metric);
    end
    baseline = baseline(keep);
    treatment = treatment(keep);

    positive = [baseline(isfinite(baseline) & baseline > 0); ...
        treatment(isfinite(treatment) & treatment > 0)];
    if isempty(positive)
        floorValue = 1e-3;
    else
        floorValue = max(min(positive) * 0.5, 1e-3);
    end
    bPlot = baseline;
    tPlot = treatment;
    bPlot(~isfinite(bPlot) | bPlot <= 0) = floorValue;
    tPlot(~isfinite(tPlot) | tPlot <= 0) = floorValue;
    bAxis = log10(bPlot);
    tAxis = log10(tPlot);

    axisRange = max([bAxis; tAxis]) - min([bAxis; tAxis]);
    pad = max(axisRange * 0.08, 0.12);
    lo = min([bAxis; tAxis]) - pad;
    hi = max([bAxis; tAxis]) + pad;
    grid = linspace(lo, hi, 220);

    fig = create_panel_figure(style.ratePanelWidthCm, style.ratePanelHeightCm);
    ax = axes(fig);
    hold(ax, 'on');
    draw_half_violin(ax, bAxis, 1, 'left', colors.baseline, grid, style);
    draw_half_violin(ax, tAxis, 2, 'right', colors.treatment, grid, style);

    rng(22 + metric_seed_offset(metric), 'twister');
    x1 = 1 + (rand(size(bAxis)) - 0.5) * 0.08;
    x2 = 2 + (rand(size(tAxis)) - 0.5) * 0.08;

    draw_categorical_change_lines(ax, x1, x2, baseline, treatment, bAxis, tAxis, colors, style);

    scatter(ax, x1, bAxis, style.channelMarkerSize, colors.baseline, 'filled', 'MarkerFaceAlpha', style.baselineDotAlpha, 'MarkerEdgeAlpha', 0);
    scatter(ax, x2, tAxis, style.channelMarkerSize, colors.treatment, 'filled', 'MarkerFaceAlpha', style.treatmentDotAlpha, 'MarkerEdgeAlpha', 0);
    draw_distribution_markers(ax, 1, baseline, floorValue, style);
    draw_distribution_markers(ax, 2, treatment, floorValue, style);

    xlim(ax, [0.52 2.48]);
    ylim(ax, [lo hi]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Baseline', treatmentLabel});
    ylabel(ax, yLabel);
    style_axis(ax, style);
    save_and_close(fig, outFile, style);
end

function offset = metric_seed_offset(metric)
    switch metric
        case 'mfr'
            offset = 0;
        case 'mbr'
            offset = 31;
        case 'ibi'
            offset = 93;
        otherwise
            offset = 0;
    end
end

function draw_categorical_change_lines(ax, x1, x2, baselineRaw, treatmentRaw, y1, y2, colors, style)
    delta = treatmentRaw(:) - baselineRaw(:);
    finite = isfinite(y1(:)) & isfinite(y2(:));
    if ~any(finite)
        return;
    end

    [~, order] = sort(abs(y2(:) - y1(:)), 'ascend', 'MissingPlacement', 'first');
    for kk = order(:)'
        if ~finite(kk)
            continue;
        end

        bActive = isfinite(baselineRaw(kk)) && baselineRaw(kk) > 0;
        tActive = isfinite(treatmentRaw(kk)) && treatmentRaw(kk) > 0;
        if ~bActive && tActive
            lineColor = colors.gained;
        elseif bActive && tActive && delta(kk) > 0
            lineColor = colors.increase;
        elseif bActive && tActive && delta(kk) < 0
            lineColor = colors.decrease;
        else
            lineColor = colors.neutralLine;
        end

        plot(ax, [x1(kk) x2(kk)], [y1(kk) y2(kk)], '-', ...
            'Color', lineColor, 'LineWidth', style.pairedDensityLineWidth);
    end
end

function draw_half_violin(ax, values, xCenter, side, color, grid, style)
    density = smooth_density(values(isfinite(values)), grid);
    if isempty(density) || max(density) <= 0
        return;
    end
    width = 0.27 * density ./ max(density);
    y = grid;
    if strcmp(side, 'left')
        fill(ax, [xCenter - width, fliplr(repmat(xCenter, size(width)))], ...
            [y, fliplr(y)], color, 'FaceAlpha', style.densityFaceAlpha, 'EdgeColor', 'none');
        plot(ax, xCenter - width, y, '-', 'Color', alpha_blend(color, style.densityLineAlpha), 'LineWidth', style.densityLineWidth);
    else
        fill(ax, [repmat(xCenter, size(width)), fliplr(xCenter + width)], ...
            [y, fliplr(y)], color, 'FaceAlpha', style.densityFaceAlpha, 'EdgeColor', 'none');
        plot(ax, xCenter + width, y, '-', 'Color', alpha_blend(color, style.densityLineAlpha), 'LineWidth', style.densityLineWidth);
    end
end

function blended = alpha_blend(color, alphaValue)
    blended = (1 - alphaValue) * [1 1 1] + alphaValue * color;
end

function density = smooth_density(values, grid)
    if numel(values) < 2
        density = zeros(size(grid));
        return;
    end
    edges = linspace(grid(1), grid(end), 37);
    counts = histcounts(values(:), edges, 'Normalization', 'pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    kernelX = linspace(-2.5, 2.5, 13);
    kernel = exp(-0.5 * kernelX .^ 2);
    kernel = kernel ./ sum(kernel);
    smooth = conv(counts, kernel, 'same');
    density = interp1(centers, smooth, grid, 'linear', 0);
end

function draw_distribution_markers(ax, x, values, floorValue, style)
    valid = values(isfinite(values));
    if isempty(valid)
        valid = floorValue;
    end
    medianValue = median(valid, 'omitnan');
    meanValue = mean(valid, 'omitnan');
    summary_marker(ax, x, medianValue, floorValue, 'median', style);
    summary_marker(ax, x, meanValue, floorValue, 'mean', style);
end

function summary_marker(ax, x, value, floorValue, kind, style)
    if isempty(value) || ~isfinite(value)
        value = floorValue;
    end
    value = log10(max(value, floorValue));
    switch kind
        case 'median'
            plot(ax, [x - 0.22, x + 0.22], [value, value], 'k-', ...
                'LineWidth', style.medianLineWidth);
        case 'mean'
            plot(ax, x, value, 'd', ...
                'MarkerSize', style.meanMarkerSize / 3, ...
                'MarkerFaceColor', 'w', ...
                'MarkerEdgeColor', 'k', ...
                'LineWidth', style.medianLineWidth);
        otherwise
            error('summary_marker:UnknownKind', 'Unknown marker kind: %s', kind);
    end
end

function plot_bursting_electrodes(data, treatmentLabel, colors, style, outFile)
    baseline = data.burstingElectrodesBaseline;
    treatment = data.burstingElectrodesTreatment;
    fig = create_panel_figure(style.ratePanelWidthCm, style.ratePanelHeightCm);
    ax = axes(fig);
    hold(ax, 'on');

    rng(211, 'twister');
    x1 = ones(size(baseline)) + (rand(size(baseline)) - 0.5) * 0.07;
    x2 = 2 * ones(size(treatment)) + (rand(size(treatment)) - 0.5) * 0.07;
    for k = 1:numel(baseline)
        if treatment(k) >= baseline(k)
            lineColor = colors.increase;
        else
            lineColor = colors.decrease;
        end
        plot(ax, [x1(k), x2(k)], [baseline(k), treatment(k)], '-', ...
            'Color', lineColor, 'LineWidth', style.pairedLineWidth);
    end
    scatter(ax, x1, baseline, style.wellMarkerSize, colors.baseline, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', style.markerEdgeWidth);
    scatter(ax, x2, treatment, style.wellMarkerSize, colors.treatment, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', style.markerEdgeWidth);
    plot(ax, [0.85, 1.15], median(baseline, 'omitnan') * [1 1], 'k-', 'LineWidth', style.medianLineWidth);
    plot(ax, [1.85, 2.15], median(treatment, 'omitnan') * [1 1], 'k-', 'LineWidth', style.medianLineWidth);

    xlim(ax, [0.55 2.45]);
    ylim(ax, [0 64]);
    set(ax, 'XTick', [1 2], 'XTickLabel', {'Baseline', treatmentLabel});
    ylabel(ax, 'Electrodes with MBR > 0');
    style_axis(ax, style);
    save_and_close(fig, outFile, style);
end

function stats = summarize_rate_features(data, study)
    stats = struct( ...
        'study', study, ...
        'metric_family', 'rate_features', ...
        'mfr_median_pct_change', median(data.mfrPct, 'omitnan'), ...
        'mbr_median_pct_change', median(data.mbrPct, 'omitnan'), ...
        'mfr_n_wells_positive', sum(data.mfrPct > 0), ...
        'mbr_n_wells_positive', sum(data.mbrPct > 0), ...
        'ibi_baseline_median_s', median(data.ibiBaseline, 'omitnan'), ...
        'ibi_treatment_median_s', median(data.ibiTreatment, 'omitnan'), ...
        'active_electrodes_baseline_median', median(data.activeElectrodesBaseline, 'omitnan'), ...
        'active_electrodes_treatment_median', median(data.activeElectrodesTreatment, 'omitnan'), ...
        'bursting_electrodes_baseline_median', median(data.burstingElectrodesBaseline, 'omitnan'), ...
        'bursting_electrodes_treatment_median', median(data.burstingElectrodesTreatment, 'omitnan'));
end

function limits = padded_limits(values, fallback)
    values = values(isfinite(values));
    if isempty(values)
        limits = fallback;
        return;
    end
    lo = min(values);
    hi = max(values);
    pad = max(10, 0.08 * (hi - lo));
    limits = [min(lo - pad, fallback(1)), max(hi + pad, 20)];
end

function style_axis(ax, style)
    box(ax, 'off');
    set(ax, ...
        'FontName', style.fontName, ...
        'FontSize', style.tickFontSize, ...
        'LineWidth', style.axesLineWidth, ...
        'TickDir', 'out', ...
        'TickLength', style.tickLength, ...
        'XGrid', 'off', ...
        'YGrid', 'off');
    ax.XAxis.FontWeight = 'normal';
    ax.YAxis.FontWeight = 'normal';
    ax.XLabel.FontName = style.fontName;
    ax.YLabel.FontName = style.fontName;
    ax.XLabel.FontSize = style.labelFontSize;
    ax.YLabel.FontSize = style.labelFontSize;
end

function save_and_close(fig, outFile, style)
    set(fig, 'Color', 'w', 'Renderer', 'painters');
    textObjects = findall(fig, '-property', 'FontName');
    for k = 1:numel(textObjects)
        try
            set(textObjects(k), 'FontName', style.fontName);
        catch
        end
    end
    save_figure(fig, outFile);
    close(fig);
    fprintf('  saved: %s\n', outFile);
end
