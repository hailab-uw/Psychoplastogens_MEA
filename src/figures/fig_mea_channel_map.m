function fig_mea_channel_map(varargin)
%FIG_MEA_CHANNEL_MAP  Publication-quality MEA grid with dual numbering.
%
%   FIG_MEA_CHANNEL_MAP() draws the 8x8 MEA grid showing:
%     - MCS electrode labels (bold, center of each cell)
%     - PZ5 channel indices (small, top-left corner)
%     - MCS 1-dim pin numbers (small, bottom-right corner)
%     - Four absent corners grayed out
%     - Internal reference electrode highlighted (MCS label 15)
%     - Color-coded by PZ5 bank (A/B/C/D)
%
%   This serves as a validation figure for the channel mapping and can
%   be included in supplementary materials or documentation.
%
%   FIG_MEA_CHANNEL_MAP('outDir', path) overrides the output directory.
%
% See also: MEA60_LAYOUT, PROJECT_CONFIG.

    cfg = project_config();

    p = inputParser;
    addParameter(p, 'outDir', output_path(cfg, '', 'shared', ''));
    parse(p, varargin{:});
    opt = p.Results;

    layout = mea60_layout();

    % Bank colors (pastel)
    bankColors = struct( ...
        'A', [0.80 0.90 1.00], ...   % light blue
        'B', [0.85 1.00 0.85], ...   % light green
        'C', [1.00 0.92 0.80], ...   % light orange
        'D', [1.00 0.85 1.00]);      % light purple
    absentColor = [0.96 0.96 0.96];  % very light gray
    refColor    = [1.00 0.95 0.80];  % light gold

    % Build lookup: for each grid position [row, col], store info
    % PZ5 channel -> bank
    function bank = get_bank(pz5ch)
        if pz5ch <= 16,      bank = 'A';
        elseif pz5ch <= 32,  bank = 'B';
        elseif pz5ch <= 48,  bank = 'C';
        else                 bank = 'D';
        end
    end

    fig = create_panel_figure(14.0, 14.0);
    ax = axes(fig, 'Position', [0.08 0.06 0.84 0.84]);
    hold(ax, 'on');

    cellW = 1;
    cellH = 1;

    % Absent corners
    absentRC = [1 1; 1 8; 8 1; 8 8];

    % Draw all 64 grid positions
    for r = 1:8
        for c = 1:8
            x = (c - 1) * cellW;
            y = (8 - r) * cellH;  % flip so row 1 is at top

            % Check if absent corner
            isAbsent = any(absentRC(:,1) == r & absentRC(:,2) == c);

            if isAbsent
                rectangle(ax, 'Position', [x y cellW cellH], ...
                    'FaceColor', absentColor, 'EdgeColor', [0.8 0.8 0.8]);
                text(ax, x + cellW/2, y + cellH/2, '--', ...
                    'HorizontalAlignment', 'center', 'FontSize', 10, ...
                    'Color', [0.7 0.7 0.7]);
                continue;
            end

            % Find which electrode is at this grid position
            mcsLabel = c * 10 + r;  % MCS convention: col*10 + row
            idx = find(layout.mcsLabels == mcsLabel);

            if isempty(idx)
                % Should not happen for non-absent positions
                rectangle(ax, 'Position', [x y cellW cellH], ...
                    'FaceColor', absentColor, 'EdgeColor', [0.8 0.8 0.8]);
                continue;
            end

            pz5ch = layout.channelList(idx);
            pin   = layout.mcsPins(idx);
            bank  = get_bank(pz5ch);

            % Choose cell color
            if mcsLabel == 15
                faceCol = refColor;
            else
                faceCol = bankColors.(bank);
            end

            % Draw cell
            rectangle(ax, 'Position', [x y cellW cellH], ...
                'FaceColor', faceCol, 'EdgeColor', [0.3 0.3 0.3], ...
                'LineWidth', 0.8);

            % MCS label (bold, center)
            labelStr = sprintf('%d', mcsLabel);
            if mcsLabel == 15
                labelStr = sprintf('%d\n(REF)', mcsLabel);
            end
            text(ax, x + cellW/2, y + cellH/2, labelStr, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 11, 'FontWeight', 'bold');

            % PZ5 channel (top-left, small)
            text(ax, x + 0.08, y + cellH - 0.12, sprintf('PZ5:%d', pz5ch), ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
                'FontSize', 5.5, 'Color', [0.3 0.3 0.6]);

            % MCS pin (bottom-right, small)
            text(ax, x + cellW - 0.08, y + 0.12, sprintf('pin %d', pin), ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
                'FontSize', 5.5, 'Color', [0.4 0.4 0.4]);
        end
    end

    % Axis labels (Arial for JNE consistency)
    for c = 1:8
        text(ax, (c - 0.5) * cellW, 8.3 * cellH, sprintf('Col %d', c), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'FontWeight', 'bold', 'FontName', 'Arial');
    end
    for r = 1:8
        text(ax, -0.4, (8 - r + 0.5) * cellH, sprintf('Row %d', r), ...
            'HorizontalAlignment', 'center', 'FontSize', 8, ...
            'FontWeight', 'bold', 'FontName', 'Arial');
    end

    % Legend for bank colors
    lgdX = 8.5;
    lgdY = 3.5;
    banks = {'A', 'B', 'C', 'D'};
    bankLabels = {'Bank A (PZ5 1-15)', 'Bank B (PZ5 17-31)', ...
                  'Bank C (PZ5 33-47)', 'Bank D (PZ5 49-63)'};
    for bi = 1:4
        rectangle(ax, 'Position', [lgdX, lgdY - (bi-1)*0.7, 0.5, 0.5], ...
            'FaceColor', bankColors.(banks{bi}), 'EdgeColor', [0.3 0.3 0.3]);
        text(ax, lgdX + 0.65, lgdY - (bi-1)*0.7 + 0.25, bankLabels{bi}, ...
            'FontSize', 6, 'VerticalAlignment', 'middle');
    end
    % REF legend
    rectangle(ax, 'Position', [lgdX, lgdY - 4*0.7, 0.5, 0.5], ...
        'FaceColor', refColor, 'EdgeColor', [0.3 0.3 0.3]);
    text(ax, lgdX + 0.65, lgdY - 4*0.7 + 0.25, 'iR Reference (pin 15)', ...
        'FontSize', 6, 'VerticalAlignment', 'middle');

    xlim(ax, [-0.8 11.5]);
    ylim(ax, [-0.2 8.6]);
    axis(ax, 'equal', 'off');

    % Save (PDF + PNG; save_figure writes both by default)
    if ~exist(opt.outDir, 'dir')
        mkdir(opt.outDir);
    end
    outFile = fullfile(opt.outDir, 'mea_channel_map');
    [pdfPath, pngPath] = save_figure(fig, outFile);
    close(fig);

    fprintf('fig_mea_channel_map: saved %s\n', pdfPath);
    fprintf('fig_mea_channel_map: saved %s\n', pngPath);
end
