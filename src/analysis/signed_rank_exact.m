function out = signed_rank_exact(baseline, treatment, varargin)
%SIGNED_RANK_EXACT Exact paired Wilcoxon signed-rank wrapper.
%
%   out = SIGNED_RANK_EXACT(baseline, treatment) returns a struct with the
%   Wilcoxon p-value, signed-rank statistic, rank-biserial effect size, and
%   a note if MATLAB cannot compute the requested exact test.

    p = inputParser;
    addRequired(p, 'baseline');
    addRequired(p, 'treatment');
    addParameter(p, 'preferExact', true, @(x) islogical(x) && isscalar(x));
    parse(p, baseline, treatment, varargin{:});
    opt = p.Results;

    b = baseline(:);
    t = treatment(:);
    if numel(b) ~= numel(t)
        error('signed_rank_exact:Args', ...
            'baseline and treatment must have equal length.');
    end

    ok = isfinite(b) & isfinite(t);
    b = b(ok);
    t = t(ok);

    out = struct('W', NaN, 'p', NaN, 'rRankBiserial', NaN, 'note', '');
    if isempty(b)
        out.note = 'no valid paired observations';
        return;
    end
    if all(t == b)
        out.W = 0;
        out.p = 1;
        out.rRankBiserial = 0;
        return;
    end
    if exist('signrank', 'file') ~= 2
        out.note = 'Statistics Toolbox (signrank) not installed';
        return;
    end

    try
        if opt.preferExact
            [wP, ~, wStats] = signrank(t, b, 'method', 'exact');
        else
            [wP, ~, wStats] = signrank(t, b);
        end
    catch ME
        try
            [wP, ~, wStats] = signrank(t, b);
            out.note = sprintf('exact signrank unavailable; used default method: %s', ME.message);
        catch ME2
            out.note = sprintf('signrank failed: %s', ME2.message);
            return;
        end
    end

    out.p = wP;
    if isfield(wStats, 'signedrank')
        out.W = wStats.signedrank;
    end
    if exist('tiedrank', 'file') == 2
        out.rRankBiserial = rank_biserial(t - b);
    elseif isempty(out.note)
        out.note = 'tiedrank unavailable; rank-biserial skipped';
    end
end

function r = rank_biserial(diffs)
    d = diffs(isfinite(diffs) & diffs ~= 0);
    if isempty(d)
        r = NaN;
        return;
    end
    ranks = tiedrank(abs(d));
    wPlus = sum(ranks(d > 0));
    wMinus = sum(ranks(d < 0));
    denom = wPlus + wMinus;
    if denom == 0
        r = 0;
    else
        r = (wPlus - wMinus) / denom;
    end
end
