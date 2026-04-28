function W = wronskian(varargin)
    % Wronskian of n vector-valued functions of t.
    %   wronskian(X)         X is an n x n sym matrix; each COLUMN is a solution x_i(t)
    %   wronskian(x1,x2,...) each xi is an n x 1 sym column
    %
    % Prints W and a verdict on linear independence.
    % Returns the simplified Wronskian if an output is requested.

    if numel(varargin) == 1
        X = varargin{1};
    else
        X = horzcat(varargin{:});
    end

    if ~isa(X, 'sym'), X = sym(X); end
    Wlocal = simplify(det(X));

    fprintf('Wronskian W(t):\n');
    disp(Wlocal);

    try
        if isAlways(Wlocal == 0, 'Unknown', 'false')
            fprintf('W = 0  =>  test inconclusive (does NOT imply dependence).\n');
        elseif isAlways(Wlocal ~= 0, 'Unknown', 'false')
            fprintf('W != 0 for all t  =>  set is linearly independent.\n');
        else
            fprintf('W depends on t. Independent at any t where W(t) != 0.\n');
        end
    catch
        fprintf('(Could not auto-classify; inspect W manually.)\n');
    end

    if nargout >= 1, W = Wlocal; end
end
