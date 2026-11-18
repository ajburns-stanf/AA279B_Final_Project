function fig = plot_grid_maps(res, prm, fig_dir, picks)
% Plot lowest-perigee maps per branch over Sun angle and halo phase
% Inputs: res - struct from grid_search, prm - struct from bcr4bp_params, fig_dir - folder for the PDF (optional, none saves nothing), picks - rows [Sun angle (deg), halo phase tau, branch p] marked T1, T2, ... on the maps (optional)
% Outputs: fig - figure handle
    if nargin < 3, fig_dir = ''; end
    if nargin < 4, picks = zeros(0, 3); end
    h0  = prm.hLEO_km;
    lab = {'Exterior', 'Interior'};
    th  = res.th_deg;  tau = res.taus;

    % Perigee maps (Parker & Anderson Figs. 3-48, 3-49): black = reaches the LEO altitude (prm.hLEO_km) or the Earth
    fig = figure('Units', 'inches', 'Position', [1 1 12 9]);  col = plot_style(fig);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    cp = [col.blue; col.vermillion; col.green; col.orange; col.purple; col.skyblue];  % T1, T2, ... colors as in plot_transfers
    mark = @(a, x, y, n) plot(a, x, y, 'o', 'MarkerSize', 12, 'MarkerFaceColor', cp(n, :), 'MarkerEdgeColor', 'w', 'LineWidth', 2.5);
    nb = 20;  % colormap rows below h0, so the black band shows on the colorbar
    for b = 1:2
        ax = nexttile(tl);  hold(ax, 'on');
        H  = log10(max(res.hmin_km(:, :, b), 1)).';
        lo = log10(h0);  hi = log10(5e5);  dl = (hi - lo)/255;
        imagesc(ax, th, tau, H, 'AlphaData', isfinite(H));          % white: no perigee before escape or Moon impact
        colormap(ax, [repmat(col.black, nb, 1); flipud(parula(255))]);  clim(ax, [lo - nb*dl, hi]);  % at or below h0 is black
        cb = colorbar(ax);  cb.Ticks = [lo 3:5];  cb.TickLabels = {num2str(h0), '$10^3$', '$10^4$', '$10^5$'};
        cb.Label.String = 'Lowest perigee alt. (km)';  cb.Label.Interpreter = 'latex';
        axis(ax, 'tight');  ax.YDir = 'normal';  ax.Layer = 'top';  xticks(ax, 0:45:360);
        title(ax, sprintf('%s branch', lab{b}));
        xlabel(ax, 'Sun angle at arrival $\theta_{s,0}$ (deg)');  ylabel(ax, 'Halo phase $\tau$');
        for n = find(picks(:, 3) == 3 - 2*b).', mark(ax, picks(n, 1), picks(n, 2), n); end
    end
    if ~isempty(picks)  % one legend for both maps, from invisible copies of the markers
        hp = arrayfun(@(n) mark(ax, NaN, NaN, n), 1:size(picks, 1));
        set(hp, {'DisplayName'}, arrayfun(@(n) sprintf('T%d', n), (1:size(picks, 1)).', 'UniformOutput', false));
        lg = legend(ax, hp, 'Orientation', 'horizontal');  lg.Layout.Tile = 'south';
    end

    if ~isempty(fig_dir)
        if ~isfolder(fig_dir), mkdir(fig_dir); end
        exportgraphics(fig, fullfile(fig_dir, 'grid_perigee_maps.pdf'), 'ContentType', 'vector');
    end
end
