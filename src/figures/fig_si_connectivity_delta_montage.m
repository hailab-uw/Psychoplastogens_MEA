function fig_si_connectivity_delta_montage(study, varargin)
%FIG_SI_CONNECTIVITY_DELTA_MONTAGE Build SI montage of per-pair delta networks.
%
%   FIG_SI_CONNECTIVITY_DELTA_MONTAGE(study) combines the per-pair
%   network_delta_labeled PNGs produced by fig_connectivity_exemplar into
%   the compact S5/S7 montage referenced by the supplementary information.

    cfg = project_config();

    p = inputParser;
    addRequired(p, 'study', @(s) any(strcmpi(s, {'doi','ket'})));
    addParameter(p, 'outDir', output_path(cfg, '', 'si', ''));
    addParameter(p, 'outName', '');
    addParameter(p, 'tileWidth', 1800);
    addParameter(p, 'pad', 60);
    parse(p, study, varargin{:});
    opt = p.Results;
    study = lower(opt.study);

    [pairs, ~] = get_pairs_and_labels(cfg, study);
    nPairs = numel(pairs);
    if nPairs == 0
        error('fig_si_connectivity_delta_montage:NoPairs', ...
            'No pairs for study %s.', study);
    end

    srcDir = fullfile(output_path(cfg, study, 'connectivity', ''), 'network');
    if ~exist(srcDir, 'dir')
        error('fig_si_connectivity_delta_montage:MissingSourceDir', ...
            'Missing connectivity network directory: %s', srcDir);
    end

    tiles = cell(nPairs, 1);
    tileHeights = zeros(nPairs, 1);
    for pi = 1:nPairs
        src = fullfile(srcDir, sprintf('%s_pair%d_network_delta_labeled.png', study, pi));
        if ~isfile(src)
            error('fig_si_connectivity_delta_montage:MissingSourceImage', ...
                'Missing source image: %s', src);
        end
        img = imread(src);
        scale = opt.tileWidth / size(img, 2);
        tile = imresize(img, scale);
        tiles{pi} = tile;
        tileHeights(pi) = size(tile, 1);
    end

    nCols = min(nPairs, 3);
    nRows = ceil(nPairs / nCols);
    tileHeight = max(tileHeights);
    pad = opt.pad;

    canvasH = nRows * tileHeight + (nRows + 1) * pad;
    canvasW = nCols * opt.tileWidth + (nCols + 1) * pad;
    canvas = uint8(255 * ones(canvasH, canvasW, 3));

    for pi = 1:nPairs
        row = floor((pi - 1) / nCols);
        col = mod(pi - 1, nCols);
        y0 = pad + row * (tileHeight + pad) + 1;
        x0 = pad + col * (opt.tileWidth + pad) + 1;

        tile = tiles{pi};
        h = size(tile, 1);
        w = size(tile, 2);
        y = y0 + floor((tileHeight - h) / 2);
        canvas(y:(y + h - 1), x0:(x0 + w - 1), :) = tile;
    end

    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    if isempty(opt.outName)
        outBase = sprintf('%s_si_connectivity_delta_networks', study);
    else
        outBase = opt.outName;
    end
    outFile = fullfile(opt.outDir, [outBase '.png']);
    imwrite(canvas, outFile);

    fprintf('fig_si_connectivity_delta_montage(%s): saved %s\n', study, outBase);
end
