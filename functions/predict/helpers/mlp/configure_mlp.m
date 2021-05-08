function layers = configure_mlp(n_inputs, n_outputs, params)
%configure_mlp
%
% Syntax: layers = configure_mlp(n_inputs, n_outputs, params)
%
% Creates an mlp as specified in params.
% Each hidden layer uses Relu activation functions
% Dropout and batch-normalisation are applied on the input and hidden
% weights.
narginchk(3, 3)
nargoutchk(1, 1)

mlp = params.mlp;

% Input layer
layers = imageInputLayer([n_inputs, 1, 1], 'Normalization', mlp.input_transform);
% Hidden layers
for layer_idx = 1:length(mlp.layer_sizes)
    % Determine sizes
    output_size = mlp.layer_sizes(layer_idx);
    
    if mlp.use_dropout
        if layer_idx == 1 % input to hidden
            dropout_value = mlp.input_dropout;
        else % hidden to hidden
            dropout_value = mlp.hidden_dropout;
        end
        layers = [layers
            dropoutLayer(dropout_value)
        ];
    end
    layers = [layers
        my_fullyConnectedLayer(output_size, ...
            'BiasLearnRateFactor', 0,...
            'WeightsInitializer', mlp.initializer)
    ];
    if mlp.use_batch_normalisation
        layers = [layers
            batchNormalizationLayer
        ];
    end
    layers = [layers
        reluLayer
    ];
end

% Output layer
layers = [layers
    my_fullyConnectedLayer(n_outputs,...
        'BiasLearnRateFactor', 0,...
        'WeightsInitializer', mlp.initializer)
];
switch mlp.error_metric
    case 'mean_squared_error'
        layers = [layers; my_meanSquaredErrorLayer()];
    case 'mean_absolute_error'
        layers = [layers; my_meanAbsoluteErrorLayer()];   
    case 'median_absolute_error'
        layers = [layers; my_medianAbsoluteErrorLayer()];   
end

end

