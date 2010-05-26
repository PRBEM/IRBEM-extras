function result = isoctave()
% result = isoctave()
% returns true in Gnu Octave
% returns false in Matlab

if exist('OCTAVE_VERSION','builtin'),
    result = true;
else
    result = false;
end
