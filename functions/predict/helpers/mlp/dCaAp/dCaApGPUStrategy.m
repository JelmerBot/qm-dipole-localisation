classdef dCaApGPUStrategy < nnet.internal.cnn.layer.util.ExecutionStrategy
    % dCaApGPUStrategy   Execution strategy for running dCaAp on the gpu
        
    methods
        function [Z, memory] = forward(~, X)
            Z = dCaAp_gpu(X);
            memory = [];
        end
        
        function [dX,dW] = backward(~, Z, dZ, X)
            dX = dCaAp_backward_gpu(Z, dZ, X);
            dW = [];
        end
    end
end