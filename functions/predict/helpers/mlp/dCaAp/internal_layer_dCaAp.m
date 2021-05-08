classdef internal_layer_dCaAp < nnet.internal.cnn.layer.Layer
    % internal_layer_dCaAp
    
    %   Copyright 2015-2017 The MathWorks, Inc.
    
    properties
        % LearnableParameters   Learnable parameters for the layer
        %   This layer has no learnable parameters.
        LearnableParameters = nnet.internal.cnn.layer.learnable.PredictionLearnableParameter.empty();
        
        % Name (char array)   A name for the layer
        Name
    end
    
    properties (Constant)
        % DefaultName   Default layer's name.
        DefaultName = 'dCaAp'
    end
    
    properties (SetAccess = private)
        % InputNames   This layer has a single input
        InputNames = {'in'}
        
        % OutputNames   This layer has a single output
        OutputNames = {'out'}
        
        % HasSizeDetermined   True for layers with size determined.
        HasSizeDetermined = true
    end
    
    properties(Access = private)
        ExecutionStrategy
    end
    
    methods
        function this = internal_layer_dCaAp(name)
            % ReLU  Constructor for the layer
            this.Name = name;
            
            this.ExecutionStrategy = dCaApHostStrategy();
        end
        
        function Z = predict( this, X )
            % predict   Forward input data through the layer and output the result
            Z = this.ExecutionStrategy.forward(X);
        end
        
        function [dX,dW] = backward( this, X, Z, dZ, ~ )
            % backward    Back propagate the derivative of the loss function
            % through the layer
            [dX,dW] = this.ExecutionStrategy.backward(Z, dZ, X);
        end
        
        function outputSize = forwardPropagateSize(~, inputSize)
            % forwardPropagateSize  Output the size of the layer based on
            %                       the input size
            outputSize = inputSize;
        end
        
        function this = inferSize(this, ~)
            % inferSize     no-op since this layer has nothing that can be
            %               inferred
        end
        
        function tf = isValidInputSize(~, ~)
            % isValidInputSize   Check if the layer can accept an input of
            % a certain size
            tf = true;
        end
        
        function this = initializeLearnableParameters(this, ~)
            % initializeLearnableParameters     no-op since there are no
            %                                   learnable parameters
        end
        
        function this = prepareForTraining(this)
            this.LearnableParameters = nnet.internal.cnn.layer.learnable.TrainingLearnableParameter.empty();
        end
        
        function this = prepareForPrediction(this)
            this.LearnableParameters = nnet.internal.cnn.layer.learnable.PredictionLearnableParameter.empty();
        end
        
        function this = setupForHostPrediction(this)
            this.ExecutionStrategy = dCaApHostStrategy();
        end
        
        function this = setupForGPUPrediction(this)
            this.ExecutionStrategy = dCaApGPUStrategy();
        end
        
        function this = setupForHostTraining(this)
            this.ExecutionStrategy = dCaApHostStrategy();
        end
        
        function this = setupForGPUTraining(this)
            this.ExecutionStrategy = dCaApGPUStrategy();
        end
    end
end
