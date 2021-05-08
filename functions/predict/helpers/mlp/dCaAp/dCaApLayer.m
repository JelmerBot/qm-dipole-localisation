function layer = dCaApLayer(varargin)
% dCaApLayer   dCaAp layer
%
%   layer = dCaApLayer() creates dCaAp layer.

% Parse the input arguments.
inputArguments = iParseInputArguments(varargin{:});

% Create an internal representation of a ReLU layer.
internalLayer = internal_layer_dCaAp(inputArguments.Name);

% Pass the internal layer to a  function to construct a user visible
% layer.
layer = layer_dCaApLayer(internalLayer);

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
inputArguments.Name = char(p.Results.Name); % make sure strings get converted to char vectors
end
