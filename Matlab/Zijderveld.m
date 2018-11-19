function [h_hor, h_ver] = Zijderveld(M, T, color_hor, color_ver)
% Plots a Zijderveld diagram for the magnetization M. 
%
% M - Nx2 or Nx3 matrix where each row is the magnetization vector at a
% given demagnetization step. 
% T - N-element vector with the demagnetization steps. Can be temperature
% but can also be AF fields. Is just used for the labels 
%
% OPTIONAL: 
% color_hor, color_ver - colors for horizontal and vertical component,
% respectively (3-element vectors)
% 
% OUTPUT: 
% h_hor - handle to the horizontal plot
% h_ver - handle to the vertical plot


    if nargin < 3
        color_hor = [0 0 0]; 
    end
    if nargin < 4
        color_ver = [0.5 0.5 0.5]; 
    end
    ho = ishold(); 
    h_hor = plot(M(:,2), M(:,1), 'o-', 'Color', color_hor);
    h_hor.MarkerFaceColor = 'w'; 
    h_ver = [];
    hold on 
    if (length(M(1,:)) == 3)
        h_ver = plot(M(:,2), M(:,3), 'o-', 'Color', color_ver);
    end
    
    if ~ho
        a = gca; 
        a.XTick = [];
        a.YTick = []; 
        a.XColor = 'w';
        a.YColor = 'w';
        grid off
        left = min(M(:,2));
        right = max(M(:,2)); 
        if size(M,2) < 3
            top = max(M(:,1));
            bottom = min(M(:,1)); 
        else
            top = max(max(M(:,1)), max(M(:,3))); 
            bottom = min(min(M(:,1)), min(M(:,3))); 
        end
        mymargin = 0.1; 
        if top<-mymargin*bottom
            top=-mymargin*bottom;
        end
        if bottom>-mymargin*top
            bottom=-mymargin*top;
        end
        if left>-mymargin*right
            left=-mymargin*right;
        end
        if right<-mymargin*left
            right=-mymargin*left;
        end
        height = top - bottom; 
        width = right - left; 
        midX = (right+left) / 2; 
        midY = (top+bottom) / 2; 
        siz = max(height, width)*1.4; 
        left   = midX - siz/2; 
        right  = midX + siz/2;
        bottom = midY - siz/2;
        top    = midY + siz/2;
        mag = 1.2;
        left2   = midX - mag*siz/2; 
        right2  = midX + mag*siz/2;
        bottom2 = midY - mag*siz/2;
        top2   = midY + mag*siz/2;
        axis([left2, right2, bottom2, top2]);
        axis('equal');
        
        siz = abs(right-left);
        ste = siz*0.02;
        plot([left right], [0 0], 'k-');
        plot([0 0], [bottom top], 'k-');
        
        if (length(M(1,:)) == 3)
            text(0, top+ste, sprintf('North (black)\n / Up (grey)'), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom'); 
            text(0, bottom-ste, sprintf('South (black)\n / Down (grey)'), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top'); 
        else
            text(0, top+ste, sprintf('N'), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom', 'FontSize', 11); 
            text(0, bottom-ste, sprintf('S'), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top', 'FontSize', 11); 
        end
        text(left-ste, 0, 'W', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Middle', 'FontSize', 11);
        text(right+ste, 0, 'E', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', 'FontSize', 11); 
    end
        
    if ho
        hold on
    else 
        hold off
    end
end