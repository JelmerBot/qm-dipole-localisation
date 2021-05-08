function c = redblue2(m)
    %REDBLUE2    Shades of red and blue color map
    %   REDBLUE2(M), is an M-by-3 matrix that defines a colormap.
    %   The colors begin with bright blue, range through shades of
    %   blue to black, and then through shades of red to bright red.
    %   REDBLUE, by itself, is the same length as the current figure's
    %   colormap. If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(redblue)
    %
    %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
    %   COLORMAP, RGBPLOT.
    %   Adam Auton, 9th October 2009
    %   Jelmer Bot, 2019
    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        b = [(0:m1-1)'/max(m1-1,1); ones(m1,1)];
        b = 1 - b;
        r = flipud(b);
        g = zeros(size(r));
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        b = [(0:m1-1)'/max(m1,1); ones(m1+1, 1)];
        b = 1 - b;
        r = flipud(b);
        g = zeros(size(r));
    end
    c = [r g b]; 
    