function velocity = normalise_input(velocity)
%normalise_input
%
% Syntax: velocity = normalise_input(velocity)
narginchk(1, 1)
nargoutchk(1, 1)
validatevelocity(velocity)

velocity.input = velocity.input ./ max(abs(velocity.input), [], 3);

end