function pre_process_figure()

set(groot, 'DefaultFigureRenderer', 'painters')

set(groot,'DefaultTextInterpreter','Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')

set(groot,'DefaultTextFontSize', 8)
set(groot, 'DefaultAxesFontSize', 8)
set(groot, 'DefaultAxesTickDir', 'out' )

set(groot, 'DefaultAxesColorOrder', ...
    [ 55 125 185;
     230 111   0;
      75 175  75;
     150  80 165;
     165  85  40;
     222  24  28] ./ 255 ...
 );

end

% Default plotting colors
% Name    Dark        Light
% red     #a01214 1   #e6191e 0.2
% blue    #275881 1   #377db9 0.2
% green   #367a34 1   #4baf4b 0.2
% purple  #6a3772 1   #9650a5 0.2
% orange  #b35900 1   #E66F00 0.2
% brown   #743c1c 1   #a55528 0.2
