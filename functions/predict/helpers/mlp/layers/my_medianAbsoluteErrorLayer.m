function layer = my_medianAbsoluteErrorLayer( varargin )
    % my_medianAbsoluteErrorLayer   Regression output layer for a neural network
    %
    %   layer = my_medianAbsoluteErrorLayer() creates a regression output layer for
    %   a neural network. The regression output layer holds the name of the
    %   loss function that is used for training the network.
    %
    %   layer = my_medianAbsoluteErrorLayer('PARAM1', VAL1) specifies optional
    %   parameter name/value pairs for creating the layer:
    %
    %       'Name'                    - A name for the layer. The default is
    %                                   ''.
    %
    %   Example:
    %       Create a regression output layer.
    %
    %       layer = my_medianAbsoluteErrorLayer();
    %
    %   See also nnet.cnn.layer.RegressionOutputLayer, classificationLayer.
    
    %   Copyright 2016-2017 The MathWorks, Inc.
    
    % Parse the input arguments
    args = iParseInputArguments(varargin{:});
    
    internalLayer = my_MedianAbsoluteError( ...
        args.Name);
    
    % Pass the internal layer to a function to construct
    layer = nnet.cnn.layer.RegressionOutputLayer(internalLayer);
    
end
    
function inputArguments = iParseInputArguments(varargin)
parser = iCreateParser();
parser.parse(varargin{:});
inputArguments = iConvertToCanonicalForm(parser);
end

function p = iCreateParser()
p = inputParser;
defaultName = '';
addParameter(p, 'Name', defaultName, @iAssertValidLayerName);
end

function iAssertValidLayerName(name)
nnet.internal.cnn.layer.paramvalidation.validateLayerName(name);
end

function inputArguments = iConvertToCanonicalForm(p)
inputArguments = struct;
inputArguments.OutputSize = [];
inputArguments.Name = char(p.Results.Name); % make sure strings get converted to char vectors
end
