function [lat_deg, lon_deg] = plotGroundTrackFromRotatingState(r_rot, vernalEq, varargin)
% plotGroundTrackFromRotatingState
%
% Input:
%   r_rot    [N x 3] or [3 x N] rotating-frame/ECEF-like position history
%   vernalEq scalar longitude offset [rad]
%
% This function:
%   1) converts rotating Cartesian position to latitude/longitude,
%   2) applies vernal equinox longitude offset,
%   3) breaks longitude discontinuities,
%   4) plots the ground track.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parse optional inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p = inputParser;

    addRequired(p,'r_rot',@(x) isnumeric(x) && ismatrix(x));
    addRequired(p,'vernalEq',@(x) isnumeric(x) && isscalar(x));

    addParameter(p,'LineWidth',1.2,@isnumeric);
    addParameter(p,'LineStyle','-',@(x) ischar(x) || isstring(x));
    addParameter(p,'Color',[],@(x) isempty(x) || (isnumeric(x) && numel(x)==3));
    addParameter(p,'NewFigure',true,@islogical);
    addParameter(p,'PlotStartStop',false,@islogical);
    addParameter(p,'DisplayName','',@(x) ischar(x) || isstring(x));
    addParameter(p,'Title','Satellite Ground Track',@(x) ischar(x) || isstring(x));

    parse(p,r_rot,vernalEq,varargin{:});

    line_width = p.Results.LineWidth;
    line_style = p.Results.LineStyle;
    line_color = p.Results.Color;
    make_new_figure = p.Results.NewFigure;
    plot_start_stop = p.Results.PlotStartStop;
    display_name = p.Results.DisplayName;
    plot_title = p.Results.Title;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ensure shape is N x 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r_rot = p.Results.r_rot;

    if size(r_rot,2) ~= 3 && size(r_rot,1) == 3
        r_rot = r_rot.';
    end

    if size(r_rot,2) ~= 3
        error('r_rot must be [N x 3] or [3 x N].');
    end

    x = r_rot(:,1);
    y = r_rot(:,2);
    z = r_rot(:,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rotating Cartesian to lat/lon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = sqrt(x.^2 + y.^2 + z.^2);

    if any(r <= 0)
        error('All position vectors must have nonzero radius.');
    end

    lat_raw_deg = asind(z ./ r);

    lon_raw_deg = rad2deg(atan2(y,x));

    % Compensate location of vernal equinox, same as your example
    lon_raw_deg = lon_raw_deg + rad2deg(vernalEq);

    % Wrap to [-180,180]
    lon_raw_deg = mod(lon_raw_deg + 180,360) - 180;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Break longitude discontinuities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lat_deg = lat_raw_deg;
    lon_deg = lon_raw_deg;

    discontinuities = find(abs(diff(lon_deg)) > 180);

    lon_deg(discontinuities + 1) = NaN;
    lat_deg(discontinuities + 1) = NaN;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if make_new_figure
        figure;
        hold on;
    else
        hold on;
    end

    % Plot coastline only once per axes
    ax = gca;
    has_coast = isappdata(ax,'GroundTrackCoastlinePlotted');

    if ~has_coast
        try
            load('topo.mat','topo');
            topoplot = [topo(:,181:360), topo(:,1:180)];
            contour(-180:179, -90:89, topoplot, [0,0], 'black', ...
                'HandleVisibility','off');
        catch
            warning('Could not load topo.mat. Plotting ground track without coastlines.');
        end

        setappdata(ax,'GroundTrackCoastlinePlotted',true);
    end

    if isempty(line_color)
        plot(lon_deg, lat_deg, ...
            'LineStyle', line_style, ...
            'LineWidth', line_width, ...
            'DisplayName', display_name);
    else
        plot(lon_deg, lat_deg, ...
            'LineStyle', line_style, ...
            'LineWidth', line_width, ...
            'Color', line_color, ...
            'DisplayName', display_name);
    end

    if plot_start_stop
        plot(lon_raw_deg(1), lat_raw_deg(1), 'bo', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor','b', ...
            'HandleVisibility','off');

        plot(lon_raw_deg(end), lat_raw_deg(end), 'rx', ...
            'MarkerSize', 8, ...
            'LineWidth', 1.5, ...
            'HandleVisibility','off');
    end

    axis equal;
    grid on;
    box on;

    xlim([-180,180]);
    xticks(-180:30:180);

    ylim([-90,90]);
    yticks(-90:30:90);

    xlabel('Longitude [deg]', 'Interpreter','latex', 'FontSize',14);
    ylabel('Latitude [deg]', 'Interpreter','latex', 'FontSize',14);
    title(plot_title, 'Interpreter','latex', 'FontSize',16);

end