function img_out = add_scale_bar(img_in, microns_per_pixel, bar_length_microns, varargin)
% ADD_SCALE_BAR  Burns a scale bar + label directly into image pixel data.
%
%   img_out = add_scale_bar(img_in, microns_per_pixel, bar_length_microns)
%   img_out = add_scale_bar(..., 'Location', 'lower-right', 'Color', [255 255 255], ...)
%
% INPUTS
%   img_in              - RGB (or grayscale) image, uint8
%   microns_per_pixel   - your calibration, e.g. 1.0 for the merged-overlay
%                          output from build_smart_warp
%   bar_length_microns  - desired physical length of the bar, e.g. 500
%
% NAME-VALUE OPTIONS (pass as 'Name', value pairs, any order)
%   'Location'    - 'lower-right' (default) | 'lower-left' | 'upper-right' | 'upper-left'
%   'Color'       - [R G B], default [255 255 255] (white)
%   'BarHeightPx' - thickness of the bar in pixels, default = 0.5% of image height
%   'Margin'      - margin from edge in pixels, default = 2% of image width
%   'ShowLabel'   - true/false, default true (draws e.g. "500 um" if insertText is available)
%   'FontScale'   - relative label size, default 1

    % ---- defaults ----
    location    = 'lower-right';
    bar_color   = [255 255 255];
    bar_height_px = -1;   % -1 = auto
    margin      = -1;     % -1 = auto
    show_label  = true;
    font_scale  = 1;

    % ---- parse name-value pairs manually (no 'arguments' block, works on old MATLAB) ----
    if mod(numel(varargin), 2) ~= 0
        error('add_scale_bar:badInput', 'Name-value arguments must come in pairs.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        switch lower(name)
            case 'location',    location = val;
            case 'color',       bar_color = val;
            case 'barheightpx', bar_height_px = val;
            case 'margin',      margin = val;
            case 'showlabel',   show_label = val;
            case 'fontscale',   font_scale = val;
            otherwise
                error('add_scale_bar:badInput', 'Unknown option "%s".', name);
        end
    end

    img_out = img_in;
    if size(img_out,3) == 1
        img_out = repmat(img_out, 1, 1, 3);   % promote grayscale to RGB so color works
    end

    [h, w, ~] = size(img_out);

    bar_length_px = round(bar_length_microns / microns_per_pixel);

    if bar_height_px < 0
        bar_height_px = max(2, round(h * 0.005));
    end

    if margin < 0
        margin = round(w * 0.02);
    end

    % compute top-left corner of the bar based on location
    switch lower(location)
        case 'lower-right'
            x0 = w - margin - bar_length_px;
            y0 = h - margin - bar_height_px;
        case 'lower-left'
            x0 = margin;
            y0 = h - margin - bar_height_px;
        case 'upper-right'
            x0 = w - margin - bar_length_px;
            y0 = margin;
        case 'upper-left'
            x0 = margin;
            y0 = margin;
        otherwise
            error('add_scale_bar:badInput', ...
                'Unknown Location "%s". Use lower-right, lower-left, upper-right, or upper-left.', location);
    end

    if x0 < 1 || y0 < 1 || (x0+bar_length_px) > w || (y0+bar_height_px) > h
        error('add_scale_bar:doesNotFit', ...
            'Scale bar does not fit within image bounds -- check microns_per_pixel and bar_length_microns.');
    end

    % draw the bar
    for c = 1:3
        img_out(y0:y0+bar_height_px-1, x0:x0+bar_length_px-1, c) = bar_color(c);
    end

    % add text label, if requested and the function is available
    if show_label
        label_str = sprintf('%g um', bar_length_microns);
        if exist('insertText', 'file') == 2 || exist('insertText', 'builtin') == 5
            font_size = max(8, round(bar_height_px * 3 * font_scale));
            text_pos_y = y0 - font_size - 4;
            if text_pos_y < 1
                text_pos_y = y0 + bar_height_px + 4;  % put below bar if no room above
            end
            img_out = insertText(img_out, [x0, text_pos_y], label_str, ...
                'FontSize', font_size, 'TextColor', bar_color, ...
                'BoxOpacity', 0, 'AnchorPoint', 'LeftTop');
        else
            warning('add_scale_bar:noTextRenderer', ...
                'insertText not available (needs Computer Vision Toolbox) -- bar drawn without text label.');
        end
    end
end
