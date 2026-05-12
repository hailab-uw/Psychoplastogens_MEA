function style = figure_style_config()
%FIGURE_STYLE_CONFIG Nature/NPP artwork parameters for generated panels.
%
%   Centralized numeric style values used directly by figure scripts.
%   Keep these explicit so each panel is clean at export time.

    style.fontName = 'Arial';
    style.tickFontSize = 6;
    style.labelFontSize = 7;
    style.titleFontSize = 7;
    style.legendFontSize = 6;
    style.panelLabelFontSize = 8;

    style.axesLineWidth = 0.6;
    style.tickLength = [0.015 0.025];
    style.referenceLineWidth = 0.45;
    style.traceLineWidth = 0.3;
    style.rasterTickWidth = 0.35;
    style.pairedLineWidth = 0.80;
    style.pairedDensityLineWidth = 0.35;
    style.densityLineWidth = 0.90;
    style.medianLineWidth = 1.10;
    style.meanMarkerSize = 18;
    style.markerEdgeWidth = 0.5;
    style.scaleBarLineWidth = 0.6;
    style.densityFaceAlpha = 0.20;
    style.densityLineAlpha = 0.80;
    style.pairedDensityPersistentAlpha = 0.40;
    style.pairedDensityStateAlpha = 0.40;
    style.baselineDotAlpha = 0.50;
    style.treatmentDotAlpha = 0.56;

    style.ratePanelWidthCm = 5.8;
    style.ratePanelHeightCm = 4.4;

    style.summaryMarkerSize = 46;
    style.channelMarkerSize = 12;
    style.wellMarkerSize = 46;

    style.baselineColor = [0.29 0.33 0.39];
    style.doiColor = [0.72 0.25 0.18];
    style.ketColor = [0.13 0.48 0.40];
    style.burstColor = [0.28 0.35 0.72];
    style.increaseColor = [0.24 0.55 0.34];
    style.decreaseColor = [0.76 0.34 0.37];
    style.silencedColor = [0.32 0.32 0.32];
    style.gainedColor = [0.24 0.39 0.78];
    style.neutralLineColor = [0.70 0.70 0.70];
    style.zeroLineColor = [0.72 0.72 0.72];
end
