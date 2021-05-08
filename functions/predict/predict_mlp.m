function predictions = predict_mlp(velocity, mlp, params)
%predict_mlp
%
% Syntax: predictions = predict_mlp(velocity, mlp, params)
%
% Implements MLP predictions.
narginchk(3, 3)
nargoutchk(1, 1)
velocity.time = 0; % force that time is removed
validatevelocity(velocity)

% Apply velocity normalisation
if params.mlp.input_normalisation
   velocity = normalise_input(velocity);
end

% Reformat for prediction    
input = permute(velocity.input, [3 2 4 1]);

% Compute predictions
predictions = predict(mlp, input);

end