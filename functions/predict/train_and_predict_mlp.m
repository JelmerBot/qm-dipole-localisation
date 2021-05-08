function [predictions, info] = train_and_predict_mlp(testing_velocity, training_velocity, params)
%train_and_predict_mlp
%
% Syntax: [predictions, info] = train_and_predict_mlp(testing_velocity, training_velocity, params)
%
% Implements MLP predictions.
narginchk(3, 3)
nargoutchk(1, 2)
testing_velocity.time = 0; % force that time is removed
training_velocity.time = 0; % force that time is removed
validatevelocity(testing_velocity)
validatevelocity(training_velocity)

% TODO
% predictions are already computed in train_mlp...
% (only when the validation set is the same as the testing set though)
% (which is always the case in this function)
[mlp, info] = train_mlp(testing_velocity, training_velocity, params);
predictions = predict_mlp(testing_velocity, mlp, params);

end