function T = export_rate_results_table(varargin)
%EXPORT_RATE_RESULTS_TABLE Manuscript summary table for rate metrics.
%
%   T = EXPORT_RATE_RESULTS_TABLE() writes a manuscript-ready table
%   for MFR, MBR, and inter-burst interval. Descriptive
%   columns follow Brofiga et al. 2023 (Q1, Q3, median, mean, SEM), while
%   inference is adapted to this paired design by testing per-well medians.

    cfg = project_config();
    p = inputParser;
    addParameter(p, 'studies', {'doi', 'ket'}, @(x) iscell(x) || isstring(x));
    addParameter(p, 'outDir', fullfile(cfg.paths.output, 'tables'), @ischar);
    addParameter(p, 'paperDir', fullfile(cfg.paths.root, 'paper', 'tables'), @ischar);
    parse(p, varargin{:});
    opt = p.Results;

    studies = cellstr(opt.studies);
    if ~exist(opt.outDir, 'dir'); mkdir(opt.outDir); end
    if ~exist(opt.paperDir, 'dir'); mkdir(opt.paperDir); end

    rows = {};
    for s = 1:numel(studies)
        study = lower(studies{s});
        [pairs, labels] = get_pairs_and_labels(cfg, study);
        specs = metric_specs();
        for mi = 1:numel(specs)
            m = specs(mi);
            accum = collect_metric(pairs, cfg.channels.recording, cfg, m);
            row = summarize_metric(study, labels.treatment, m, accum);
            rows{end + 1, 1} = row; %#ok<AGROW>
        end
    end

    T = struct2table(vertcat(rows{:}));
    T.p_wilcoxon_bh_rate_family = nan(height(T), 1);
    studyNames = unique(T.study);
    for i = 1:numel(studyNames)
        studyMask = T.study == studyNames(i);
        if studyNames(i) ~= "doi"
            continue;
        end
        pVals = T.p_wilcoxon_well_pct_change(studyMask);
        [pAdj, ~] = bh_fdr(pVals, cfg.stats.fdr_q);
        T.p_wilcoxon_bh_rate_family(studyMask) = pAdj(:);
    end
    csvPath = fullfile(opt.outDir, 'rate_results_table.csv');
    mdPath = fullfile(opt.outDir, 'rate_results_table.md');
    paperMdPath = fullfile(opt.paperDir, 'rate_results_table.md');
    writetable(T, csvPath);
    write_markdown_table(T, mdPath);
    write_markdown_table(T, paperMdPath);
    fprintf('Wrote %s\n', csvPath);
    fprintf('Wrote %s\n', mdPath);
    fprintf('Wrote %s\n', paperMdPath);
end

function specs = metric_specs()
    specs = [
        struct('id', 'mfr', 'label', 'Mean firing rate', ...
            'abbrev', 'MFR', 'unit', 'spikes/min')
        struct('id', 'mbr', 'label', 'Mean bursting rate', ...
            'abbrev', 'MBR', 'unit', 'bursts/min')
        struct('id', 'ibi', 'label', 'Inter-burst interval', ...
            'abbrev', 'IBI', 'unit', 's')
    ];
end

function accum = collect_metric(pairs, channels, cfg, spec)
    accum.poolB = [];
    accum.poolT = [];
    accum.wellB = nan(numel(pairs), 1);
    accum.wellT = nan(numel(pairs), 1);
    accum.wellPct = nan(numel(pairs), 1);

    for k = 1:numel(pairs)
        bCache = load_cache(pairs(k).baseline, cfg);
        tCache = load_cache(pairs(k).treatment, cfg);
        switch spec.id
            case 'mfr'
                [b, t] = align_cache_metric(bCache, tCache, channels, 'spikeRates');
                keep = b >= cfg.silent.min_rate_spike | t >= cfg.silent.min_rate_spike;
                pctKeep = b >= cfg.silent.min_rate_spike;
            case 'mbr'
                [b, t] = align_cache_metric(bCache, tCache, channels, 'burstRates');
                keep = b > 0 | t > 0;
                pctKeep = b > 0;
            case 'ibi'
                [b, t] = align_burst_feature(bCache, tCache, channels, 'ibi');
                keep = isfinite(b) & isfinite(t) & b > 0 & t > 0;
                pctKeep = keep;
            otherwise
                error('collect_metric:UnknownMetric', 'Unknown metric: %s', spec.id);
        end
        pctVals = safe_pct(b(pctKeep), t(pctKeep));
        b = b(keep);
        t = t(keep);
        accum.poolB = [accum.poolB; b(:)];
        accum.poolT = [accum.poolT; t(:)];
        accum.wellB(k) = median(b, 'omitnan');
        accum.wellT(k) = median(t, 'omitnan');
        accum.wellPct(k) = median(pctVals, 'omitnan');
    end
end

function row = summarize_metric(study, treatmentLabel, spec, accum)
    b = accum.poolB(:);
    t = accum.poolT(:);
    wellOK = isfinite(accum.wellB) & isfinite(accum.wellT);
    wellB = accum.wellB(wellOK);
    wellT = accum.wellT(wellOK);

    row = struct();
    row.study = string(study);
    row.treatment = string(treatmentLabel);
    row.metric = string(spec.abbrev);
    row.metric_label = string(spec.label);
    row.unit = string(spec.unit);
    row.n_wells = numel(wellB);
    row.n_electrode_pairs = numel(b);

    row.baseline_q1 = percentile_value(b, 25);
    row.baseline_q3 = percentile_value(b, 75);
    row.baseline_median = median(b, 'omitnan');
    row.baseline_mean = mean(b, 'omitnan');
    row.baseline_sem = sem_value(b);

    row.treatment_q1 = percentile_value(t, 25);
    row.treatment_q3 = percentile_value(t, 75);
    row.treatment_median = median(t, 'omitnan');
    row.treatment_mean = mean(t, 'omitnan');
    row.treatment_sem = sem_value(t);

    wellPct = accum.wellPct(isfinite(accum.wellPct));
    row.per_well_median_delta = median(wellT - wellB, 'omitnan');
    row.per_well_median_pct_change = median(wellPct, 'omitnan');
    row.n_wells_increased = sum(wellPct > 0);
    row.p_wilcoxon_well_pct_change = signed_rank_p(zeros(size(wellPct)), wellPct, true);
    row.p_wilcoxon_paired_electrodes_descriptive = signed_rank_p(b, t, false);
end

function [b, t] = align_cache_metric(bCache, tCache, channels, metric)
    chB = bCache.channelsUsed(:)';
    chT = tCache.channelsUsed(:)';
    [~, idxB] = ismember(channels, chB);
    [~, idxT] = ismember(channels, chT);
    ok = idxB > 0 & idxT > 0;
    bAll = bCache.(metric)(:);
    tAll = tCache.(metric)(:);
    b = bAll(idxB(ok));
    t = tAll(idxT(ok));
end

function [b, t] = align_burst_feature(bCache, tCache, channels, feature)
    chB = bCache.channelsUsed(:)';
    chT = tCache.channelsUsed(:)';
    [~, idxB] = ismember(channels, chB);
    [~, idxT] = ismember(channels, chT);
    ok = idxB > 0 & idxT > 0;
    idxB = idxB(ok);
    idxT = idxT(ok);
    b = nan(numel(idxB), 1);
    t = nan(numel(idxT), 1);
    for i = 1:numel(idxB)
        b(i) = channel_burst_feature(bCache, idxB(i), feature);
        t(i) = channel_burst_feature(tCache, idxT(i), feature);
    end
end

function value = channel_burst_feature(cacheStruct, idx, feature)
    starts = cacheStruct.burstStartTimes{idx}(:);
    ends = cacheStruct.burstEndTimes{idx}(:);
    switch feature
        case 'duration'
            n = min(numel(starts), numel(ends));
            if n == 0
                value = NaN;
                return;
            end
            values = 1000 * (ends(1:n) - starts(1:n));
        case 'ibi'
            starts = starts(isfinite(starts));
            if numel(starts) < 2
                value = NaN;
                return;
            end
            values = diff(sort(starts));
        otherwise
            error('channel_burst_feature:UnknownFeature', ...
                'Unknown feature: %s', feature);
    end
    values = values(isfinite(values) & values > 0);
    value = median(values, 'omitnan');
end

function value = percentile_value(x, p)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = prctile(x, p);
    end
end

function value = sem_value(x)
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = std(x, 0) / sqrt(numel(x));
    end
end

function p = signed_rank_p(b, t, preferExact)
    stat = signed_rank_exact(b, t, 'preferExact', preferExact);
    p = stat.p;
end

function pct = safe_pct(b, t)
    pct = nan(size(b));
    nz = b ~= 0;
    pct(nz) = 100 * (t(nz) - b(nz)) ./ b(nz);
    pct = pct(isfinite(pct));
end

function write_markdown_table(T, path)
    fid = fopen(path, 'w');
    if fid < 0
        error('write_markdown_table:OpenFailed', 'Could not open %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '| Study | Metric | Unit | n wells | n electrode pairs | Baseline Q1 | Baseline Q3 | Baseline median | Baseline mean | Baseline SEM | Treatment Q1 | Treatment Q3 | Treatment median | Treatment mean | Treatment SEM | Median delta | Median %% change | Wilcoxon p, well %% change | BH p, DOI rate family | Wilcoxon p, paired electrodes |\n');
    fprintf(fid, '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
    for i = 1:height(T)
        fprintf(fid, '| %s | %s | %s | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
            char(T.treatment(i)), char(T.metric(i)), char(T.unit(i)), ...
            T.n_wells(i), T.n_electrode_pairs(i), ...
            fmt(T.baseline_q1(i)), fmt(T.baseline_q3(i)), ...
            fmt(T.baseline_median(i)), fmt(T.baseline_mean(i)), fmt(T.baseline_sem(i)), ...
            fmt(T.treatment_q1(i)), fmt(T.treatment_q3(i)), ...
            fmt(T.treatment_median(i)), fmt(T.treatment_mean(i)), fmt(T.treatment_sem(i)), ...
            fmt(T.per_well_median_delta(i)), fmt(T.per_well_median_pct_change(i)), ...
            fmt_p(T.p_wilcoxon_well_pct_change(i)), ...
            fmt_p(T.p_wilcoxon_bh_rate_family(i)), ...
            fmt_p(T.p_wilcoxon_paired_electrodes_descriptive(i)));
    end
    fprintf(fid, '\nDescriptive columns are pooled electrode-level values using the same inclusion rules as the rate-feature panels. The primary inferential column is the paired Wilcoxon signed-rank test on per-well median percent change; BH correction is applied only across the DOI rate family (MFR, MBR, IBI). The paired-electrode p-value is included only as a descriptive comparison to the plotted electrode-level distributions.\n');
end

function s = fmt(x)
    if ~isfinite(x)
        s = 'NA';
    elseif abs(x) >= 100
        s = sprintf('%.1f', x);
    elseif abs(x) >= 10
        s = sprintf('%.2f', x);
    else
        s = sprintf('%.3f', x);
    end
end

function s = fmt_p(x)
    if ~isfinite(x)
        s = 'NA';
    elseif x < 0.001
        s = '<0.001';
    else
        s = sprintf('%.3f', x);
    end
end
