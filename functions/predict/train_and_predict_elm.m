function [predictions, train_error, train_predictions] = train_and_predict_elm(testing_velocity, training_velocity, params)
%train_and_predict_elm
%
% Syntax: [predictions, train_error, train_predictions] = train_and_predict_elm(testing_velocity, training_velocity, params)
%
% Implements ELM predictions.
narginchk(3, 3)
nargoutchk(1, 3)
testing_velocity.time = 0;
validatevelocity(testing_velocity)
training_velocity.time = 0;
validatevelocity(training_velocity)
if size(testing_velocity.input, 3) ~= size(training_velocity.input, 3)
    error('Testing and training inputs have to have the same size per source!')
end

elm = train_elm(training_velocity, params);
predictions = predict_elm(testing_velocity, elm, params);            

% Also compute training error if requested
if nargout >= 2
    train_predictions = predict_elm(training_velocity, elm, params);
    train_error = mean(compute_prediction_errors(train_predictions,...
                                                 training_velocity.sources),...
                                                 params);
end
end

