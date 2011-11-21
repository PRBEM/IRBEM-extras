function varargout = kdtree(what,varargin)
% tree = kdtree('build',X) % build kdtree
%
% [index,R2,XNN] = kdtree('kNN',X,tree,X0,k,DistScale) % retrieve index and (optionally location XNN) of k points nearest to X0
% DistScale (size(X0))is optional - defines weighting for each dimension
%
% [index,R2,XNN] = kdtree('multikNN',X,tree,X0,k,DistScale) % retrieve index and (optionally location XNN) of k points nearest to X0
% ... = kdtree('multikNN',...,'max') % use maximum squared difference
%           rather than sum squared difference for distance measure
%           (default is 'sum')
% DistScale (length size(X0,2))is optional - defines weighting for each dimension
% optionally handles multiple X0, where each is a row vector
% size(index) and size(R2) are [k size(X0,1)]
% size(XNN) is [k size(X0,1) size(X0,2)]
% requires kdtree.dll (or .so) and kdtree.h in the matlab path
% R2 is squared distance

varargout = cell(1,nargout);

switch(lower(what)),
    case {'nolibbuild'},
        [varargout{:}] = kdtree_build(varargin{1});
    case {'libbuild','build'},
        [varargout{:}] = kdtree_libbuild(varargin{1});
    case {'knn'},
        [varargout{:}] = kdtree_kNN(varargin{:});
    case {'multiknn'},
        [varargout{:}] = kdtree_multikNN(varargin{:});
    otherwise
        error('Unknown what "%s"',what);
end

end

function tree = kdtree_build(X)

% size variables
[Nx,Nc] = size(X);
if (Nx==0),
    error('Cannot create empty kdtree');
end

% save memory, use smallest necessary int type for index records
itype = get_int_type(Nx);
ctype = get_int_type(Nc);

tree.root = itype(0);
tree.c = ctype(zeros(Nx,1));
tree.parent = itype(zeros(Nx,1));
tree.left = itype(zeros(Nx,1));
tree.right = itype(zeros(Nx,1));

[foo,tree] = kdtree_split(X,tree,itype(1:size(X,1))',ctype(1),itype(0));

end % end of function


function tree = kdtree_libbuild(X)
% size variables
[Nx,Nc] = size(X);
if Nx==0,
    error('Cannot create empty kdtree');
end
flags = 0; % Matlab is column major (0)
kdtree_loadlib;

if isoctave,
    [tree.root,tree.c,tree.parent,tree.left,tree.right] = kdtree_build_oct(X,Nx,Nc,flags);
else
    
    % int kdtree_build(const double *X,
    % 		 const unsigned long int Nx,
    % 		 const unsigned long int Nc,
    % 		 const int flags,
    % 		 unsigned long int *root,
    % 		 unsigned short int *c,
    % 		 unsigned long int *parent,
    % 		 unsigned long int *left,
    % 		 unsigned long int *right);
    
    
    
    
    rootPtr = libpointer('ulongPtr',0);
    cPtr = libpointer('uint16Ptr',zeros(Nx,1)); % assumes ushort = uint16, should be true
    parentPtr = libpointer('ulongPtr',zeros(Nx,1));
    leftPtr = libpointer('ulongPtr',zeros(Nx,1));
    rightPtr = libpointer('ulongPtr',zeros(Nx,1));
    
    calllib('kdtree','kdtree_build',X,Nx,Nc,flags,rootPtr,cPtr,parentPtr,leftPtr,rightPtr);
    
    tree.root = rootPtr.value;
    tree.c = cPtr.value;
    tree.parent = parentPtr.value;
    tree.left = leftPtr.value;
    tree.right = rightPtr.value;
end

end % end of function

function [leaf_index,tree] = kdtree_split(X,tree,I,c,iparent)

[Nx,Nc] = size(X);
NI = size(I,1);

[xc,isort] = sort(X(I,c));
I = I(isort);

jmid = ceil(NI/2);
leaf_index = I(jmid);

if iparent==0,
    tree.root = leaf_index;
end

tree.c(leaf_index) = c;
tree.parent(leaf_index) = iparent;

% prepare to call split
c = 1+rem(c,Nc);
i1 = (1:(jmid-1));
i2 = ((jmid+1):NI);

if ~isempty(i1),
    [leaf,tree] = kdtree_split(X,tree,I(i1),c,leaf_index);
    tree.left(leaf_index) = leaf;
end
if ~isempty(i2),
    [leaf,tree] = kdtree_split(X,tree,I(i2),c,leaf_index);
    tree.right(leaf_index) = leaf;
end
end % end of function


function varargout = kdtree_kNN(X,tree,X0,k,DistScale)

if nargin < 5,
    DistScale = ones(1,size(X0,2));
    dist_fun = @(x,y) sum((x(:)-y(:)).^2);
else
    dist_fun = @(x,y) sum((DistScale(:).*(x(:)-y(:))).^2);
end

varargout = cell(1,nargout);

persistent warned
if isoctave,
    if isempty(warned) || ~warned,
        warning('***kdtree_kNN does not work in octave because it uses nested functions. Using multikNN instead');
        warned = true;
    end
    [varargout{:}] = kdtree_multikNN(X,tree,X0,k,DistScale);
    return
end


% variables shared with subfunctions: X, tree, Nc, index, XNN, R2, dist_fun, DistScale
Nc = size(X0,2);
index = zeros(k,1,class(tree.parent(1)));
R2 = inf(k,1);
if nargout >= 3,
    XNN = nan(k,Nc);
else
    XNN = [];
end

% first we descend to closest point
i = tree.root;
while i>0,
    ilast = i;
    if X0(tree.c(i))>X(i,tree.c(i)),
        i = tree.right(i);
    else
        i = tree.left(i);
    end
end

% now ilast points to the closest leaf to X0

% check parent nodes for closer children
check_node(ilast,true,true,true);

varargout{1} = index;
if nargout >= 2,
    varargout{2} = R2;
end
if nargout >= 3,
    varargout{3} = XNN;
end

    function check_node(i,check_parent,check_left,check_right)
        if i <= 0, % null node
            return;
        end
        dxc = DistScale(tree.c(i))*(X(i,tree.c(i))-X0(tree.c(i)));
        dxc2 = dxc^2;
        if dxc2 <= R2(end), % this node and its children could be nearest neighbors
            r2i = dist_fun(X(i,:),X0);
            if r2i <= R2(end), % this is an active nearest neighbor
                iR2 = length(R2);
                while (iR2 > 1) && (r2i <= R2(iR2-1)),
                    R2(iR2) = R2(iR2-1);
                    index(iR2) = index(iR2-1);
                    if ~isempty(XNN),
                        XNN(iR2,:) = XNN(iR2-1,:);
                    end
                    iR2 = iR2 - 1;
                end
                R2(iR2) = r2i;
                index(iR2) = i;
                if ~isempty(XNN),
                    XNN(iR2,:) = X(i,:);
                end
            end
            if check_right,
                check_node(tree.right(i),false,true,true);
            end
            if check_left,
                check_node(tree.left(i),false,true,true);
            end
        else
            % this node is too far away, but its kids might be
            if dxc<0, % only kids to the right can have smaller dxc2
                if check_right,
                    check_node(tree.right(i),false,true,true);
                end
            else % only kids to the left can have smaller dxc2
                if check_left,
                    check_node(tree.left(i),false,true,true);
                end
            end
        end
        if check_parent && (tree.parent(i)>0),
            % check parent, but don't recheck this node
            isright = (tree.right(tree.parent(i))==i);
            check_node(tree.parent(i),true,isright,~isright);
        end
    end

end

function [index,R2,XNN] = kdtree_multikNN(X,tree,X0,k,varargin)

DistScale = [];
how = 'sum'; % how to combine differences to compute distance
if length(varargin)==2,
    DistScale = varargin{1};
    how = varargin{2};
elseif length(varargin)==1,
    if ischar(varargin{1}),
        how = varargin{1};
    else
        DistScale = varargin{1};
    end
end

switch(lower(how)),
    case 'max',
        how_flag = 2;
    case 'sum',
        how_flag = 0;
    otherwise
        error('how="%s" not understood',how);
end

[Nx,Nc] = size(X);
NX0 = size(X0,1);
flags = 0; % Matlab is column major (0)
flags = flags + how_flag; % add flag for how to compbine differences into dist
kdtree_loadlib;

if isoctave,
    [index,R2] = kdtree_multikNN_oct(X,Nx,Nc,flags,tree.root,tree.c,tree.parent,tree.left,tree.right,X0,NX0,k,DistScale);
    index = index';
    R2 = R2';
else % end of octave specific code
    
    if isempty(DistScale),
        DistScale = libpointer('doublePtr'); % null
    end
    
    indexPtr = libpointer('ulongPtr',zeros(NX0,k));
    if nargout >= 2,
        R2Ptr = libpointer('doublePtr',zeros(NX0,k));
    else
        R2Ptr = libpointer('doublePtr'); % null
    end
    
    %int32 kdtree_kNN(doublePtr, ulong, ulong, int32, ulong, uint16Ptr, ulongPtr, ulongPtr, ulongPtr, doublePtr, ulong, ulong, doublePtr, ulongPtr, doublePtr)
    

    calllib('kdtree','kdtree_kNN',X,Nx,Nc,flags,...
        tree.root,tree.c,tree.parent,tree.left,tree.right,...
        X0,NX0,k,DistScale,indexPtr,R2Ptr);
    
    index = get(indexPtr,'val')';
    if nargout >= 2,
        R2 = get(R2Ptr,'val')';
    end
    
end % end of matlab-specific code

if nargout >= 3, % populate XNN output
    XNN = nan([k size(X0)]);
    for c = 1:Nc,
        XNN(:,:,c) = reshape(X(index,c),size(index));
    end
end

end % end of function


function itype = get_int_type(N)

Ibits = ceil(log(1+N)/log(2));

if Ibits <= 8,
    itype = @uint8;
elseif Ibits <= 16,
    itype = @uint16;
elseif Ibits <= 32,
    itype = @uint32;
elseif Ibits <= 64,
    itype = @uint64;
end

end

function kdtree_loadlib()
if isoctave,
    if ~exist('kdtree_multikNN_oct'), % load
        octfile = file_in_loadpath('kdtree_oct.oct');
        if isempty(octfile),
            error('unable to locate kdtree_oct.oct');
        end
        autoload('kdtree_multikNN_oct',octfile);
        autoload('kdtree_build_oct',octfile);
    end
else
    if ~libisloaded('kdtree'),
        
        if ~exist('kdtree.h','file'),
            error('Unable to locate kdtree.h');
        end
        
        if ispc,
            libname = 'kdtree.dll';
        else
            libname = 'kdtree.so';
        end
        
        loadlibrary(libname,'kdtree.h','alias','kdtree');
        
    end
end
end
