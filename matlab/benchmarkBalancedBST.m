%% benchmarkBalancedBST
%  Compare mr.aux.BalancedBST against MATLAB's dictionary() for
%  insert and lookup of string-keyed scalar data.
%
%  Requires R2022b+ for dictionary().

%% ---- parameters --------------------------------------------------------
Ns = [100, 500, 1000, 5000, 10000];%, 50000, 100000];  % number of entries to test
nTrials = 5;                                        % repeat & take median

fprintf('%-8s | %-26s | %-26s | %-26s | %-26s\n', ...
    'N', 'BST insert (ms)', 'map insert (ms)', ...
    'BST lookup (ms)', 'map lookup (ms)');
fprintf('%s\n', repmat('-', 1, 124));

for ni = 1:numel(Ns)
    N = Ns(ni);

    % ---- generate random keys & values ---------------------------------
    keys = cell(1, N);
    vals = 1:N;
    for k = 1:N
        keys{k} = sprintf('key_%16d', randi(1e8));
    end
    % ensure uniqueness (duplicates would skew insert count)
    [keys, ia] = unique(keys, 'stable');
    vals = vals(ia);
    N = numel(keys);    
    
    % look-up keys: 50 % present, 50 % absent
    nLookup = N;
    lookupKeys = cell(1, nLookup);
    for k = 1:nLookup
        if rand < 0.5
            lookupKeys{k} = keys{randi(N)};        % existing key
        else
            lookupKeys{k} = sprintf('miss_%16d', randi(1e8));
        end
    end

    % % test for strings() -- actually slower
    % for k = 1:N
    %     keys{k} = string(keys{k});
    %     lookupKeys{k} = string(lookupKeys{k});
    % end

    t_bst_ins  = zeros(1, nTrials);
    t_dict_ins = zeros(1, nTrials);
    t_bst_lkp  = zeros(1, nTrials);
    t_dict_lkp = zeros(1, nTrials);

    for trial = 1:nTrials
        % ================================================================
        %  BalancedBST insert
        % ================================================================
        tree = mr.aux.BalancedBST();
        tic;
        for k = 1:N
            insert(tree, keys{k}, vals(k)); % avoid tree. syntax to avoid overloaded subsref/subsasgn calls
        end
        t_bst_ins(trial) = toc;

        % ================================================================
        %  dictionary insert
        % ================================================================
        d = containers.Map('KeyType','char','ValueType','double');%dictionary(string.empty, double.empty);
        tic;
        for k = 1:N
            d(keys{k}) = vals(k);
        end
        t_dict_ins(trial) = toc;

        % ================================================================
        %  BalancedBST lookup
        % ================================================================
        tic;
        for k = 1:nLookup
            [~, ~] = lookup(tree, lookupKeys{k}, 0); % avoid tree. syntax to avoid overloaded subsref/subsasgn calls
        end
        t_bst_lkp(trial) = toc;

        % ================================================================
        %  dictionary lookup
        % ================================================================
        tic;
        for k = 1:nLookup
            if isKey(d, lookupKeys{k}) %isKey(d, lookupKeys{k}) is ~ 2.5 times faster than d.isKey(lookupKeys{k})
                tmp = d(lookupKeys{k}); 
            end
        end
        t_dict_lkp(trial) = toc;
    end

    fprintf('%-8d | %10.2f  (%6.2f us/op) | %10.2f  (%6.2f us/op) | %10.2f  (%6.2f us/op) | %10.2f  (%6.2f us/op)\n', ...
        N, ...
        median(t_bst_ins)*1e3,  median(t_bst_ins)/N*1e6, ...
        median(t_dict_ins)*1e3, median(t_dict_ins)/N*1e6, ...
        median(t_bst_lkp)*1e3,  median(t_bst_lkp)/nLookup*1e6, ...
        median(t_dict_lkp)*1e3, median(t_dict_lkp)/nLookup*1e6);
end

%% ---- correctness smoke test -------------------------------------------
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('Correctness smoke test ... ');
tree = mr.aux.BalancedBST();
ref  = containers.Map('KeyType','char','ValueType','double');%dictionary(string.empty, double.empty);
rng(42);
testKeys = cell(1, 500);
for k = 1:500
    testKeys{k} = sprintf('k%d', randi(200));
end
for k = 1:500
    v = randi(9999);
    tree.insert(testKeys{k}, v);
    ref(testKeys{k}) = v;
end
allKeys = ref.keys();
ok = true;
for k = 1:numel(allKeys)
    [val, index] = tree.lookup(allKeys{k}, -1);
    if 0==index || val ~= ref(allKeys{k})
        ok = false;
        fprintf('MISMATCH at key "%s"\n', allKeys{k});
    end
end
% check a missing key
[val, index] = tree.lookup('__nonexistent__', -999);
if 0~=index || val ~= -999
    ok = false;
    fprintf('MISMATCH for missing key\n');
end
if ok
    fprintf('PASSED\n');
else
    fprintf('FAILED\n');
end

%% ---- update correctness test ----------------------------------------
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('update() correctness test ... ');
tree2 = mr.aux.BalancedBST();
tree2.insert('alpha', 1);
tree2.insert('beta',  2);
tree2.insert('gamma', 3);
[~, idx] = tree2.lookup('beta', 0);
tree2.update(idx, 42);
[val, ~] = tree2.lookup('beta', 0);
ok2 = (val == 42);
% make sure other keys are intact
[v1, ~] = tree2.lookup('alpha', 0); ok2 = ok2 && (v1 == 1);
[v3, ~] = tree2.lookup('gamma', 0); ok2 = ok2 && (v3 == 3);
if ok2, fprintf('PASSED\n'); else, fprintf('FAILED\n'); end

%% ---- remove correctness test -------------------------------------------
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('remove() correctness test ... ');
tree3 = mr.aux.BalancedBST();
ref3  = containers.Map('KeyType','char','ValueType','double');
rng(99);
nIns = 1000;
tkeys = cell(1, nIns);
for k = 1:nIns
    tkeys{k} = sprintf('r%d', randi(400));
end
for k = 1:nIns
    v = randi(9999);
    tree3.insert(tkeys{k}, v);
    ref3(tkeys{k}) = v;
end
% remove ~half the keys
allK = ref3.keys();
removeList = allK(randperm(numel(allK), floor(numel(allK)/2)));
for k = 1:numel(removeList)
    r1 = tree3.remove(removeList{k});
    remove(ref3, removeList{k});
    if ~r1
        ok = false;
        fprintf('remove() returned false for existing key "%s"\n', removeList{k});
    end
end
% verify remaining keys
remainK = ref3.keys();
ok3 = true;
for k = 1:numel(remainK)
    [val, idx] = tree3.lookup(remainK{k}, -1);
    if idx == 0 || val ~= ref3(remainK{k})
        ok3 = false;
        fprintf('MISMATCH at key "%s" after remove\n', remainK{k});
    end
end
% verify removed keys are gone
for k = 1:numel(removeList)
    [~, idx] = tree3.lookup(removeList{k}, -1);
    if idx ~= 0
        ok3 = false;
        fprintf('Key "%s" still present after remove\n', removeList{k});
    end
end
% verify length
if tree3.length() ~= numel(remainK)
    ok3 = false;
    fprintf('length() mismatch: got %d, expected %d\n', tree3.length(), numel(remainK));
end
% verify remove of non-existent key returns false
if tree3.remove('__nope__')
    ok3 = false;
    fprintf('remove() returned true for non-existent key\n');
end
if ok3, fprintf('PASSED\n'); else, fprintf('FAILED\n'); end

%% ---- remove + re-insert (free list reuse) test -------------------------
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('remove + re-insert (free-list reuse) test ... ');
tree4 = mr.aux.BalancedBST();
for k = 1:100
    tree4.insert(sprintf('x%d', k), k);
end
for k = 1:50
    tree4.remove(sprintf('x%d', k));
end
% re-insert — these should reuse freed slots
for k = 1:50
    tree4.insert(sprintf('y%d', k), k + 1000);
end
ok4 = (tree4.length() == 100);
for k = 51:100
    [v, idx] = tree4.lookup(sprintf('x%d', k), -1);
    if idx == 0 || v ~= k, ok4 = false; break; end
end
for k = 1:50
    [v, idx] = tree4.lookup(sprintf('y%d', k), -1);
    if idx == 0 || v ~= k + 1000, ok4 = false; break; end
end
if ok4, fprintf('PASSED\n'); else, fprintf('FAILED\n'); end

%% ---- remove benchmark -------------------------------------------------
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('remove() benchmark\n');
fprintf('%-8s | %-26s | %-26s\n', ...
    'N', 'BST remove (ms)', 'map remove (ms)');
fprintf('%s\n', repmat('-', 1, 65));

for ni = 1:numel(Ns)
    N = Ns(ni);
    keys4 = cell(1, N);
    for k = 1:N
        keys4{k} = sprintf('key_%16d', k);  % unique keys
    end

    t_bst_rm  = zeros(1, nTrials);
    t_dict_rm = zeros(1, nTrials);

    for trial = 1:nTrials
        % build BST
        tree5 = mr.aux.BalancedBST();
        for k = 1:N, tree5.insert(keys4{k}, k); end
        % remove all
        perm = randperm(N);
        tic;
        for k = 1:N
            remove(tree5, keys4{perm(k)}); % avoid overloaded subref calls
        end
        t_bst_rm(trial) = toc;

        % build map
        d5 = containers.Map('KeyType','char','ValueType','double');
        for k = 1:N, d5(keys4{k}) = k; end
        tic;
        for k = 1:N
            remove(d5, keys4{perm(k)});
        end
        t_dict_rm(trial) = toc;
    end

    fprintf('%-8d | %10.2f  (%6.2f us/op) | %10.2f  (%6.2f us/op)\n', ...
        N, ...
        median(t_bst_rm)*1e3,  median(t_bst_rm)/N*1e6, ...
        median(t_dict_rm)*1e3, median(t_dict_rm)/N*1e6);
end
