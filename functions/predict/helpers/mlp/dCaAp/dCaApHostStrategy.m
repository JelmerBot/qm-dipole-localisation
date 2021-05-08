classdef dCaApHostStrategy < nnet.internal.cnn.layer.util.ExecutionStrategy
    % dCaApHostStrategy   Execution strategy for running dCaAp on the host
        
    methods
        function [Z, memory] = forward(~, X)
            Z = dCaAp(X);
            memory = [];
        end
        
        function [dX,dW] = backward(~, Z, dZ, X)
            dX = dCaAp_backward(Z, dZ, X);
            dW = [];
        end
    end
end