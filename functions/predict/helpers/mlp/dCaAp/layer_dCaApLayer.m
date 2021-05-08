classdef layer_dCaApLayer < nnet.cnn.layer.Layer & nnet.internal.cnn.layer.Externalizable
    % layer_dCaApLayer
    %
    % Used by dCaApLayer.
    
    properties(Dependent)
        % Name   A name for the layer
        %   The name for the layer. If this is set to '', then a name will
        %   be automatically set at training time.
        Name         
    end
    
    methods
        function this = layer_dCaApLayer(privateLayer)
            this.PrivateLayer = privateLayer;
        end      
        
        function val = get.Name(this)
            val = this.PrivateLayer.Name;
        end
        
        function this = set.Name(this, val)
            iAssertValidLayerName(val);
            this.PrivateLayer.Name = char(val);
        end
        
        function out = saveobj(this)
            out.Version = 1.0;
            out.Name = this.PrivateLayer.Name;
        end
    end
    
    methods(Static)
        function this = loadobj(in)
            internalLayer = internal_layer_dCaAp(in.Name);
            this = layer_dCaApLayer(internalLayer);
        end
    end

    methods(Access = protected)
        function [description, type] = getOneLineDisplay(~)
            description = iGetMessageString('nnet_cnn:layer:dCaApLayer:oneLineDisplay');
            
            type = iGetMessageString( 'nnet_cnn:layer:dCaApLayer:Type' );
        end
    end
end

function messageString = iGetMessageString( varargin )
messageString = getString( message( varargin{:} ) );
end

function iAssertValidLayerName(name)
nnet.internal.cnn.layer.paramvalidation.validateLayerName(name);
end