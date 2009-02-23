%% make test files for var2cdf

pathstr = fileparts(which(mfilename));

%% first test: a 1-d vector
vector = (1:10);
% vector =
%      1     2     3     4     5     6     7     8     9    10

var2cdf([pathstr,filesep,'test_vector.cdf'],vector);

%% second test: a 2-d matrix

matrix = reshape(1:12,[3 4]);
% matrix =
%      1     4     7    10
%      2     5     8    11
%      3     6     9    12

var2cdf([pathstr,filesep,'test_matrix.cdf'],matrix);

%% third test: a 3-d tensor
tensor = reshape(1:24,[2 3 4]);
% tensor(:,:,1) =
%      1     3     5
%      2     4     6
% tensor(:,:,2) =
%      7     9    11
%      8    10    12
% tensor(:,:,3) =
%     13    15    17
%     14    16    18
% tensor(:,:,4) =
% 
%     19    21    23
%     20    22    24
var2cdf([pathstr,filesep,'test_tensor.cdf'],tensor);

%% fourth test: a complicated, deeply nested structure
structure = struct('name','sample');
structure.vector = vector;
structure.matrix = matrix;
structure.tensor = tensor;
structure.substructure = structure;
structure.substructure.subsubstructure = structure;

var2cdf([pathstr,filesep,'test_structure.cdf'],structure);

% structure.name
% sample
% structure.vector
%      1     2     3     4     5     6     7     8     9    10
% structure.matrix
%      1     4     7    10
%      2     5     8    11
%      3     6     9    12
% structure.tensor
% (:,:,1) =
%      1     3     5
%      2     4     6
% (:,:,2) =
%      7     9    11
%      8    10    12
% (:,:,3) =
%     13    15    17
%     14    16    18
% (:,:,4) =
%     19    21    23
%     20    22    24
% structure.substructure.name
% sample
% structure.substructure.vector
%      1     2     3     4     5     6     7     8     9    10
% structure.substructure.matrix
%      1     4     7    10
%      2     5     8    11
%      3     6     9    12
% structure.substructure.tensor
% (:,:,1) =
%      1     3     5
%      2     4     6
% (:,:,2) =
%      7     9    11
%      8    10    12
% (:,:,3) =
%     13    15    17
%     14    16    18
% (:,:,4) =
%     19    21    23
%     20    22    24
% structure.substructure.subsubstructure.name
% sample
% structure.substructure.subsubstructure.vector
%      1     2     3     4     5     6     7     8     9    10
% structure.substructure.subsubstructure.matrix
%      1     4     7    10
%      2     5     8    11
%      3     6     9    12
% structure.substructure.subsubstructure.tensor
% (:,:,1) =
%      1     3     5
%      2     4     6
% (:,:,2) =
%      7     9    11
%      8    10    12
% (:,:,3) =
%     13    15    17
%     14    16    18
% (:,:,4) =
%     19    21    23
%     20    22    24
% structure.substructure.subsubstructure.substructure.name
% sample
% structure.substructure.subsubstructure.substructure.vector
%      1     2     3     4     5     6     7     8     9    10
% structure.substructure.subsubstructure.substructure.matrix
%      1     4     7    10
%      2     5     8    11
%      3     6     9    12
% structure.substructure.subsubstructure.substructure.tensor
% (:,:,1) =
%      1     3     5
%      2     4     6
% (:,:,2) =
%      7     9    11
%      8    10    12
% (:,:,3) =
%     13    15    17
%     14    16    18
% (:,:,4) =
%     19    21    23
%    

%% this little code snippet will print the structure
% x = {'structure',structure};
% i=1;
% while i < length(x),
%     base = x{i};
%     s = x{i+1};
%     fnames = fieldnames(s);
%     for j = 1:length(fnames),
%         var = fnames{j};
%         if isstruct(s.(var)),
%             x = {x{:},[base,'.',var],s.(var)};
%         else
%             disp(sprintf('%s.%s',base,var));
%             disp(s.(var));
%         end
%     end
%     i = i+2;
% end


