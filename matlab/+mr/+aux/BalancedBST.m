classdef BalancedBST < handle
    %BALANCEDBST  AVL-balanced binary search tree mapping keys to data. See
    %             https://en.wikipedia.org/wiki/AVL_tree for further infos.
    %
    %   Keys are typically char vectors (strings); data is typically a
    %   scalar, but both can be any MATLAB type whose comparison is
    %   well-defined (numeric or char/string for keys).
    %
    %   Performance notes
    %   -----------------
    %   Cell arrays in MATLAB are contiguous pointer arrays with O(1)
    %   random access — NOT linked lists.  The tree traversal is
    %   O(log n) key-comparisons.  The pre-allocated parallel-array
    %   storage avoids per-node heap allocations; doubling on overflow
    %   gives amortised O(1) growth.
    %
    %   Usage
    %   -----
    %     tree = mr.aux.BalancedBST();
    %     tree.insert('alpha', 1);
    %     tree.insert('beta',  2);
    %
    %     [val, index] = tree.lookup('alpha', 0);   % val=1, index=non-zero
    %     [val, index] = tree.lookup('gamma', 0);   % val=0, index=0
    %

    % --- storage: parallel arrays (pre-allocated, grown as needed) -------
    properties (Access = private)
        keys      % cell(1,cap)  – stored keys
        vals      % cell(1,cap)  – stored data values
        L         % int32(1,cap) – left-child index  (0 = none)
        R         % int32(1,cap) – right-child index (0 = none)
        H         % int32(1,cap) – subtree height
        rootIdx   % int32        – root node index   (0 = empty tree)
        cnt       % int32        – number of allocated slots (high-water mark)
        cap       % int32        – allocated capacity
        keyIsChar % bool         – whether the key is a char string or an arbitrary vector
        pathBuf   % int32
        dirBuf    % int32
        freeHead  % int32        – head of free-slot linked list (0 = none)
        liveCnt   % int32        – number of live (non-deleted) nodes
    end

    % =====================================================================
    %  PUBLIC INTERFACE
    % =====================================================================
    methods (Access = public)

        function obj = BalancedBST(initialCapacity)
            %BALANCEDBST  Construct an empty tree.
            %   tree = mr.aux.BalancedBST()
            %   tree = mr.aux.BalancedBST(initialCapacity)
            if nargin < 1, initialCapacity = 64; end
            obj.cap       = int32(initialCapacity);
            obj.keys      = cell(1, obj.cap);
            obj.vals      = cell(1, obj.cap);
            obj.L         = zeros(1, obj.cap, 'int32');
            obj.R         = zeros(1, obj.cap, 'int32');
            obj.H         = zeros(1, obj.cap, 'int32');
            obj.rootIdx   = int32(0);
            obj.cnt       = int32(0);
            obj.keyIsChar = false;
            obj.pathBuf   = int32(0);
            obj.dirBuf    = int32(0);
            obj.freeHead  = int32(0);
            obj.liveCnt   = int32(0);
        end

        function [val, index] = lookup(obj, key, default)
            %LOOKUP  Search for KEY in the tree.
            %   [val, index]       = tree.lookup(key, default)
            %
            %   If KEY is present, VAL is the associated data and
            %   INDEX is non-zero. Otherwise VAL = DEFAULT and INDEX is
            %   zero.
            %

            val   = default;
            index = 0;

            % --- cache property arrays as locals (COW = free) ------------
            %   Every obj.prop(idx) in a loop pays ~300 ns dispatch
            %   overhead; local variable access is ~1 ns.  For a
            %   read-only traversal COW means the snapshot is free.
            %kk = obj.keys;
            %LL = obj.L;
            %RR = obj.R;

            idx = obj.rootIdx;
            while idx ~= int32(0)
                c = mr.aux.BalancedBST.compareKeys(key, obj.keys{idx});
                if c == 0
                    val   = obj.vals{idx};   % single property access
                    index = idx;
                    break;
                end
                if c < 0
                    idx = obj.L(idx);
                else
                    idx = obj.R(idx);
                end
            end
        end

        function insert(obj, key, val)
            %INSERT  Insert or update a key-value pair.
            %   tree.insert(key, val)
            %
            %   If KEY already exists its data is overwritten.

            root = obj.rootIdx;
            if root == int32(0) % the tree was empty up to now
                obj.rootIdx = obj.allocNode(key, val);
                obj.keyIsChar = ischar(key);
                return;
            end

            % --- search phase --------------------------------------------
            depth = int32(0);
            idx   = root;

            while idx ~= int32(0)
                c = mr.aux.BalancedBST.compareKeys(key, obj.keys{idx});
                depth = depth + 1;
                obj.pathBuf(depth) = idx;
                if c == 0
                    obj.vals{idx} = val;     % overwrite existing value
                    return;
                elseif c < 0
                    obj.dirBuf(depth) = int32(-1);
                    idx = obj.L(idx);
                else
                    obj.dirBuf(depth) = int32(1);
                    idx = obj.R(idx);
                end
            end

            % --- allocate new node and link to parent --------------------
            newIdx = obj.allocNode(key, val);
            pIdx = obj.pathBuf(depth);
            if obj.dirBuf(depth) < 0
                obj.L(pIdx) = newIdx;
            else
                obj.R(pIdx) = newIdx;
            end

            % --- rebalance bottom-up -------------------------------------
            obj.rebalanceUp(depth);
        end

        function update(obj, index, val)
            %UPDATE  Update the value at a known node index in O(1).
            %   tree.update(index, newVal)
            %
            %   INDEX is typically obtained from a prior tree.lookup().
            %   The key is unchanged; only the associated data is
            %   overwritten.  No rebalancing is needed.
            obj.vals{index} = val;
        end

        function removed = remove(obj, key)
            %REMOVE  Remove the entry with the given key.
            %   removed = tree.remove(key)
            %
            %   Returns true if KEY was found and removed, false if
            %   KEY was not present in the tree.

            if obj.rootIdx == int32(0)
                removed = false;
                return;
            end

            % --- search for the node -------------------------------------
            depth = int32(0);
            idx   = obj.rootIdx;
            while idx ~= int32(0)
                c = mr.aux.BalancedBST.compareKeys(key, obj.keys{idx});
                depth = depth + 1;
                obj.pathBuf(depth) = idx;
                if c == 0
                    break;
                elseif c < 0
                    obj.dirBuf(depth) = int32(-1);
                    idx = obj.L(idx);
                else
                    obj.dirBuf(depth) = int32(1);
                    idx = obj.R(idx);
                end
            end

            if idx == int32(0)
                removed = false;
                return;
            end

            % idx == obj.pathBuf(depth) is the node to delete
            li = obj.L(idx);
            ri = obj.R(idx);

            if li ~= 0 && ri ~= 0
                % --- TWO CHILDREN: replace with in-order successor --------
                %   Go right once, then left as far as possible.
                obj.dirBuf(depth) = int32(1);   % going right from idx
                succ = ri;
                depth = depth + 1;
                obj.pathBuf(depth) = succ;
                while obj.L(succ) ~= int32(0)
                    obj.dirBuf(depth) = int32(-1);
                    succ = obj.L(succ);
                    depth = depth + 1;
                    obj.pathBuf(depth) = succ;
                end
                % Copy successor's key/value to the target node
                obj.keys{idx} = obj.keys{succ};
                obj.vals{idx} = obj.vals{succ};
                % Successor has at most a right child
                replacement = obj.R(succ);
                obj.freeNode(succ);
            else
                % --- ZERO or ONE CHILD -----------------------------------
                if li ~= 0
                    replacement = li;
                else
                    replacement = ri;   % may be 0 (leaf)
                end
                obj.freeNode(idx);
            end

            % --- link replacement to parent of deleted node ---------------
            if depth > 1
                p = obj.pathBuf(depth - 1);
                if obj.dirBuf(depth - 1) < 0
                    obj.L(p) = replacement;
                else
                    obj.R(p) = replacement;
                end
            else
                obj.rootIdx = replacement;
            end

            % --- rebalance from parent of deleted node upward -------------
            obj.rebalanceUp(depth - 1);

            removed = true;
        end

        function n = length(obj)
            %LENGTH  Number of live key-value pairs in the tree.
            n = double(obj.liveCnt);
        end

       % function varargout = subsref(obj, S)
       %     %SUBSREF  Overloaded subscript reference.
       %     %   val = tree('key')   — equivalent to tree.lookup('key', [])
       %     %
       %     %   Dot-reference (tree.method, tree.prop) and curly-brace
       %     %   indexing are forwarded to the built-in handler so that
       %     %   normal method calls keep working.
       %     if S(1).type(1) == '('
       %         key = S(1).subs{1};
       %         [val, ~] = obj.lookup(key, []);
       %         if numel(S) > 1
       %             % chained indexing, e.g. tree('key').field
       %             [varargout{1:nargout}] = subsref(val, S(2:end));
       %         else
       %             varargout{1} = val;
       %         end
       %     else
       %         % '.' or '{}' — delegate to built-in
       %         [varargout{1:nargout}] = builtin('subsref', obj, S);
       %     end
       % end

        function obj = subsasgn(obj, S, val)
            %SUBSASGN  Overloaded subscript assignment.
            %   tree('key') = val   — equivalent to tree.insert('key', val)
            %
            %   Dot-assignment and curly-brace assignment are forwarded
            %   to the built-in handler.
            if S(1).type(1) == '(' && numel(S) == 1
                key = S(1).subs{1};
                obj.insert(key, val);
            else
                % '.' or '{}' or chained — delegate to built-in
                obj = builtin('subsasgn', obj, S, val);
            end
        end

    end % public methods

    % =====================================================================
    %  PRIVATE HELPERS
    % =====================================================================
    methods (Access = private)

        % ----- node allocation -------------------------------------------
        function idx = allocNode(obj, key, val)
            if obj.freeHead ~= int32(0)
                idx = obj.freeHead;
                obj.freeHead = obj.L(idx);  % L was reused as next-free
            else
                obj.cnt = obj.cnt + 1;
                if obj.cnt > obj.cap
                    obj.grow();
                end
                idx = obj.cnt;
            end
            obj.keys{idx} = key;
            obj.vals{idx} = val;
            obj.L(idx)    = int32(0);
            obj.R(idx)    = int32(0);
            obj.H(idx)    = int32(1);
            obj.liveCnt   = obj.liveCnt + 1;
        end

        function freeNode(obj, idx)
            %FREENODE  Return a slot to the free list.
            obj.keys{idx} = [];
            obj.vals{idx} = [];
            obj.R(idx)    = int32(0);
            obj.H(idx)    = int32(0);
            obj.L(idx)    = obj.freeHead;   % reuse L as next-free pointer
            obj.freeHead  = idx;
            obj.liveCnt   = obj.liveCnt - 1;
        end

        function grow(obj)
            added  = obj.cap;            % double the capacity
            obj.keys = [obj.keys, cell(1, added)];
            obj.vals = [obj.vals, cell(1, added)];
            obj.L    = [obj.L,    zeros(1, added, 'int32')];
            obj.R    = [obj.R,    zeros(1, added, 'int32')];
            obj.H    = [obj.H,    zeros(1, added, 'int32')];
            obj.cap  = obj.cap + int32(added);
        end

        % ----- rebalance from pathBuf(depth) up to root ------------------
        function rebalanceUp(obj, depth)
            for i = depth:-1:1
                nd = obj.pathBuf(i);

                % -- refresh height (inlined) --
                li = obj.L(nd);  ri = obj.R(nd);
                lh = int32(0);   rh = int32(0);
                if li ~= 0, lh = obj.H(li); end
                if ri ~= 0, rh = obj.H(ri); end
                obj.H(nd) = int32(1) + max(lh, rh);

                bf  = rh - lh;
                nnd = nd;

                if bf < -1
                    % left-heavy
                    child = obj.L(nd);
                    cli = obj.L(child);  cri = obj.R(child);
                    clh = int32(0);      crh = int32(0);
                    if cli ~= 0, clh = obj.H(cli); end
                    if cri ~= 0, crh = obj.H(cri); end
                    if (crh - clh) > 0    % Left-Right case
                        gc = obj.R(child);
                        obj.R(child) = obj.L(gc);
                        obj.L(gc)    = child;
                        tl = obj.L(child); tr = obj.R(child);
                        tlh = int32(0); trh = int32(0);
                        if tl ~= 0, tlh = obj.H(tl); end
                        if tr ~= 0, trh = obj.H(tr); end
                        obj.H(child) = int32(1) + max(tlh, trh);
                        tl = obj.L(gc); tr = obj.R(gc);
                        tlh = int32(0); trh = int32(0);
                        if tl ~= 0, tlh = obj.H(tl); end
                        if tr ~= 0, trh = obj.H(tr); end
                        obj.H(gc) = int32(1) + max(tlh, trh);
                        obj.L(nd) = gc;
                    end
                    x = obj.L(nd);
                    obj.L(nd) = obj.R(x);
                    obj.R(x)  = nd;
                    tl = obj.L(nd); tr = obj.R(nd);
                    tlh = int32(0); trh = int32(0);
                    if tl ~= 0, tlh = obj.H(tl); end
                    if tr ~= 0, trh = obj.H(tr); end
                    obj.H(nd) = int32(1) + max(tlh, trh);
                    tl = obj.L(x); tr = obj.R(x);
                    tlh = int32(0); trh = int32(0);
                    if tl ~= 0, tlh = obj.H(tl); end
                    if tr ~= 0, trh = obj.H(tr); end
                    obj.H(x) = int32(1) + max(tlh, trh);
                    nnd = x;

                elseif bf > 1
                    % right-heavy
                    child = obj.R(nd);
                    cli = obj.L(child);  cri = obj.R(child);
                    clh = int32(0);      crh = int32(0);
                    if cli ~= 0, clh = obj.H(cli); end
                    if cri ~= 0, crh = obj.H(cri); end
                    if (crh - clh) < 0    % Right-Left case
                        gc = obj.L(child);
                        obj.L(child) = obj.R(gc);
                        obj.R(gc)    = child;
                        tl = obj.L(child); tr = obj.R(child);
                        tlh = int32(0); trh = int32(0);
                        if tl ~= 0, tlh = obj.H(tl); end
                        if tr ~= 0, trh = obj.H(tr); end
                        obj.H(child) = int32(1) + max(tlh, trh);
                        tl = obj.L(gc); tr = obj.R(gc);
                        tlh = int32(0); trh = int32(0);
                        if tl ~= 0, tlh = obj.H(tl); end
                        if tr ~= 0, trh = obj.H(tr); end
                        obj.H(gc) = int32(1) + max(tlh, trh);
                        obj.R(nd) = gc;
                    end
                    y = obj.R(nd);
                    obj.R(nd) = obj.L(y);
                    obj.L(y)  = nd;
                    tl = obj.L(nd); tr = obj.R(nd);
                    tlh = int32(0); trh = int32(0);
                    if tl ~= 0, tlh = obj.H(tl); end
                    if tr ~= 0, trh = obj.H(tr); end
                    obj.H(nd) = int32(1) + max(tlh, trh);
                    tl = obj.L(y); tr = obj.R(y);
                    tlh = int32(0); trh = int32(0);
                    if tl ~= 0, tlh = obj.H(tl); end
                    if tr ~= 0, trh = obj.H(tr); end
                    obj.H(y) = int32(1) + max(tlh, trh);
                    nnd = y;
                end

                % -- re-link to parent ------------------------------------
                if nnd ~= nd
                    if i > 1
                        p = obj.pathBuf(i-1);
                        if obj.dirBuf(i-1) < 0
                            obj.L(p) = nnd;
                        else
                            obj.R(p) = nnd;
                        end
                    else
                        obj.rootIdx = nnd;
                    end
                end
            end
        end

    end % private methods

    % =====================================================================
    %  STATIC (key comparison)
    % =====================================================================
    methods (Static, Access = private)

        function c = compareKeys(a, b)
            %COMPAREKEYS  Lexicographic comparison returning -1, 0, or +1.
            %   Handles numeric keys (scalar <, >, ==) and char/string
            %   keys (character-by-character comparison).
            % [~,I]=sort({a,b});
            % c = diff(I)*~strcmp(a,b);
            %la = numel(a); lb = numel(b);
            %ml = min(la, lb);
            for k = 1:min(numel(a),numel(b))
                % c=sign(int32(a(k))-int32(b(k)));
                % if c~=0
                %     return;
                % end
                if a(k) > b(k)
                    c = 1;
                    return;
                elseif a(k) < b(k)
                    c = -1;
                    return;
                end
            end
            c = sign(numel(a)-numel(b));
            % a=string(a);
            % b=string(b);
            % if a>b
            %     c = 1;
            % elseif a < b
            %     c = -1;
            % else
            %     c = 0;
            % end
        end

    end % static methods
end
