function [C, names] = movement_regressors(rparam_file, mr_types)
% returns functions of movement parameters, and regressor names 

if nargin < 1
  rparam_file = spm_get([0 1], 'rp*.txt', ...
			'Select movement parameter file');
end
if nargin < 2
  mr_types = {'moves'};
end

% Names for movement parameters
move_names = {'x trans', 'y trans', 'z trans', ...
	       'x rot', 'y rot', 'z rot'};
cmn = char(move_names);

moves = spm_load(rparam_file);
if isempty(moves)
  error(['Cannot get movement parameters from: ' rparam_file]);
end
if size(moves, 2) ~= 6
  error(['Expecting a Nx6 movement parameter matrix from file' ...
	rparam_file]);
end

% mean centre moves
moves = moves - ones(size(moves, 1), 1) * mean(moves);
moves_m1 = [zeros(1, 6); moves(1:end-1,:)];

C = []; 
names = {}; 
    
if ismember('moves', mr_types)
  C = [C moves];
  names = [names move_names];
end
if ismember('mm1', mr_types)
  C = [C moves_m1];
  names = [names sf_name_cat(move_names, ' minus 1')];
end
if ismember('moves_2', mr_types)
  C = [C moves.^2];
  names = [names sf_name_cat(move_names, ' .^2')];
end
if ismember('mm1_2', mr_types)
  C = [C moves_m1.^2];
  names = [names sf_name_cat(move_names, ' minus 1 .^2')];
end
return

function names = sf_name_cat(n, suffix)
names = cellstr([char(n) repmat(suffix, 6, 1)])';
return
