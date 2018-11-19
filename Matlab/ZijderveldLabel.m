function ZijderveldLabel(M, T, labels, color_hor, location, color_ver)
% Adds labels to a Zijderveld diagram. 
%
% M - Nx2 or Nx3 matrix where each row is the magnetization vector at a
% given demagnetization step. 
% T - N-element vector with the demagnetization steps. Can be temperature
% but can also be AF fields. Is just used for the label positioning 
% labels - add labels to the steps given in this N-element vector. 
% 
% OPTIONAL
% color_hor, color_ver - colors for horizontal and vertical component,
% respectively (3-element vectors)
% location - 'topleft' or 'bottomright', determining where the labels are 
% placed relative to the markers. Default: 'bottomright'

    if nargin < 4 
        color_hor = [0 0 0]; 
    end
    if nargin < 5
        location = 'bottomright';
    end
    if nargin < 6
        color_ver = color_hor; 
    end
    
    font_size = 10;

    ax = axis;
    wid = ax(2)-ax(1);
    if strcmpi(location, 'bottomright')
        w = wid * 0.02 *sin(-pi/4);
        h = wid * 0.02 *cos(-pi/4); 
        align = 'left'; 
    else
        w = -wid * 0.02 *sin(-pi/4);
        h = -wid * 0.02 *cos(-pi/4); 
        align = 'right'; 
    end
    lastx = Inf; 
    lasty = Inf; 
    siz = getpixelposition(gca);
    wid_pix = siz(3);
    for n = 1:length(labels)
        id = find(T>=labels(n) & ~isnan(M(:,1)') & ~isnan(M(:,2)'), 1, 'first');
        x = M(id,2)-w; 
        y = M(id,1)-h; 
        d = sqrt((x-lastx)^2 + (y-lasty)^2);
        if d/wid*wid_pix > 10
            text(x, y, strcat(' ', num2str(labels(n))), ...
                    'Rotation',-45, 'FontSize', font_size, 'Color', color_hor, ...
                    'HorizontalAlignment', align); 
            lastx = x;
            lasty = y; 
            if size(M, 2) >= 3
                text(M(id,2)-w, M(id,3)-h, strcat(' ', num2str(labels(n))), ...
                        'Rotation',-45, 'FontSize', font_size, 'Color', color_ver, ...
                        'HorizontalAlignment', align); 
            end
        end
    end
    
end