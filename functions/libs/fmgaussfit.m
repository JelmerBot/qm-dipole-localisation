function [fitresult, resnorm] = fmgaussfit(xx,yy,zz, params)
% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface 
%   height. XX and YY are vectors or matrices defining the x and y 
%   components of a surface. If XX and YY are vectors, length(XX) = n and 
%   length(YY) = m, where [m,n] = size(Z). In this case, the vertices of the
%   surface faces are (XX(j), YY(i), ZZ(i,j)) triples. To create XX and YY 
%   matrices for arbitrary domains, use the meshgrid function. FMGAUSSFIT
%   uses the lsqcurvefit tool, and the OPTIMZATION TOOLBOX. The initial
%   guess for the gaussian is places at the maxima in the ZZ plane. The fit
%   is restricted to be in the span of XX and YY.
%   See:
%       http://en.wikipedia.org/wiki/Gaussian_function
%          
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, resnorm] = fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.
%   Copyright 2013, Nathan Orloff.
%   Adapted by Jelmer Bot, 2019.
%% Condition the data
[xData, yData, zData] = prepareSurfaceData( xx, yy, zz );
xyData = {xData,yData};
%% Set up the startpoint
[amp, ind] = max(zData); % amp is the amplitude.
xo = xData(ind); % guess that it is at the maximum
yo = yData(ind); % guess that it is at the maximum
ang = 45; % angle in degrees.
sy = 1;
sx = 1;
zo = median(zData(:))-std(zData(:));
xmax = max(params.domain.x_range);
ymax = max(params.domain.y_range);
xmin = min(params.domain.x_range);
ymin = min(params.domain.y_range);
%% Set up fittype and options.
Lower = [0, 0, 0, 0, xmin, ymin, 0];
Upper = [Inf, 180, Inf, Inf, xmax, ymax, Inf]; % angles greater than 90 are redundant
StartPoint = [amp, ang, sx, sy, xo, yo, zo];%[amp, sx, sxy, sy, xo, yo, zo];
options = optimset('Algorithm','trust-region-reflective',...
    'Display','off',...
    'MaxFunEvals',2e3,...
    'MaxIter',100,...
    'TolX',1e-3,...
    'TolFun',eps,...
    'TolCon',eps);

%% perform the fitting
[fitresult,resnorm] = ...
    lsqcurvefit(@gaussian2D,StartPoint,xyData,zData,Lower,Upper,options);
end