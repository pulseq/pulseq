function benchmarkMapPerformance()
%BENCHMARKMAPPERFORMANCE Compare performance of struct, containers.Map,
%   Java HashMap, and dictionary for string-keyed storage.
%
%   Tests insertion, random lookup, overwrite, key-existence check, deletion,
%   and iteration across multiple collection sizes.
%
%   The dictionary object requires MATLAB R2022b or later. If unavailable,
%   dictionary tests are skipped automatically.

    sizes = [100, 500, 1000 5000, 10000];
    nTrials = 7;  % repeat each measurement and take the median
    keyLen  = 12; % length of each random-string key
    testJava = true;
    testDict = hasDictionary();

    nImpl = 4; % struct, Map, Java, dictionary

    headers = {'N', ...
               'Struct_Insert','Map_Insert','Java_Insert','Dict_Insert', ...
               'Struct_Lookup','Map_Lookup','Java_Lookup','Dict_Lookup', ...
               'Struct_Lookup2','Map_Lookup2','Java_Lookup2','Dict_Lookup2', ...
               'Struct_Overwrite','Map_Overwrite','Java_Overwrite','Dict_Overwrite', ...
               'Struct_HasKey','Map_HasKey','Java_HasKey','Dict_HasKey', ...
               'Struct_HasKey2','Map_HasKey2','Java_HasKey2','Dict_HasKey2', ...
               'Struct_Iterate','Map_Iterate','Java_Iterate','Dict_Iterate', ...
               'Struct_Delete','Map_Delete','Java_Delete','Dict_Delete'};

    results = zeros(numel(sizes), numel(headers));

    for si = 1:numel(sizes)
        N = sizes(si);
        fprintf('\n========== N = %d ==========\n', N);

        % ------ generate random keys (valid MATLAB identifiers) ----------
        keys = generateRandomKeys(N, keyLen);
        keys2 = generateRandomKeys(N, keyLen);
        vals = randi(1e6, 1, N);           % random integer values
        lookupOrder  = randperm(N);         % randomised access order
        overwriteVals = randi(1e6, 1, N);   % new values for overwrite test

        % Pre-allocate timing arrays
        tInsert   = zeros(nTrials, nImpl);
        tLookup   = zeros(nTrials, nImpl);
        tLookup2  = zeros(nTrials, nImpl);
        tOverwrite= zeros(nTrials, nImpl);
        tHasKey   = zeros(nTrials, nImpl);
        tHasKey2  = zeros(nTrials, nImpl);
        tIterate  = zeros(nTrials, nImpl);
        tDelete   = zeros(nTrials, nImpl);

        for t = 1:nTrials
            % ============================================================
            %  1. INSERTION
            % ============================================================

            % --- struct ---
            S = struct();
            tic;
            for i = 1:N
                S.(keys{i}) = vals(i);
            end
            tInsert(t,1) = toc;

            % --- containers.Map ---
            tic;
            %M = containers.Map(keys, vals);
            M = containers.Map();
            for i = 1:N
                M(keys{i})=vals(i);
            end
            tInsert(t,2) = toc;

            % --- java.util.HashMap ---
            if testJava
                tic;
                J = java.util.HashMap(N);
                for i = 1:N
                    J.put(keys{i}, vals(i));
                end
                tInsert(t,3) = toc;
            end

            % --- dictionary ---
            if testDict
                tic;
                %D = dictionary(string(keys), vals);
                D = dictionary();
                for i = 1:N
                    D(keys{i})=vals(i);
                end
                tInsert(t,4) = toc;
            end

            % ============================================================
            %  2. RANDOM LOOKUP EXISTING
            % ============================================================
            dummy = 0; %#ok<NASGU> % prevent optimiser from removing the loop

            % --- struct ---
            tic;
            for i = lookupOrder
                dummy = S.(keys{i});
            end
            tLookup(t,1) = toc;

            % --- containers.Map ---
            tic;
            for i = lookupOrder
                dummy = M(keys{i});
            end
            tLookup(t,2) = toc;

            % --- java.util.HashMap ---
            if testJava
                tic;
                for i = lookupOrder
                    dummy = J.get(keys{i});
                end
                tLookup(t,3) = toc;
            end

            % --- dictionary ---
            if testDict
                tic;
                for i = lookupOrder
                    dummy = D(keys{i});
                end
                tLookup(t,4) = toc;
            end

            % ============================================================
            %  3. RANDOM LOOKUP EXCEPTION
            % ============================================================
            dummy = 0; %#ok<NASGU> % prevent optimiser from removing the loop

            % --- struct ---
            tic;
            for i = lookupOrder
                try
                    dummy = S.(keys2{i});
                catch
                end
            end
            tLookup2(t,1) = toc;

            % --- containers.Map ---
            tic;
            for i = lookupOrder
                try
                    dummy = M(keys2{i});
                catch
                end
            end
            tLookup2(t,2) = toc;

            % --- java.util.HashMap ---
            if testJava
                tic;
                for i = lookupOrder
                    try
                        dummy = J.get(keys2{i});
                    catch
                    end
                end
                tLookup2(t,3) = toc;
            end

            % --- dictionary ---
            if testDict
                tic;
                for i = lookupOrder
                    try
                        dummy = D(keys2{i});
                    catch
                    end
                end
                tLookup2(t,4) = toc;
            end

            % ============================================================
            %  4. OVERWRITE (update existing keys)
            % ============================================================

            % --- struct ---
            tic;
            for i = 1:N
                S.(keys{i}) = overwriteVals(i);
            end
            tOverwrite(t,1) = toc;

            % --- containers.Map ---
            tic;
            for i = 1:N
                M(keys{i}) = overwriteVals(i);
            end
            tOverwrite(t,2) = toc;

            % --- java.util.HashMap ---
            if testJava
                tic;
                for i = 1:N
                    J.put(keys{i}, overwriteVals(i));
                end
                tOverwrite(t,3) = toc;
            end

            % --- dictionary ---
            if testDict
                tic;
                for i = 1:N
                    D(keys{i}) = overwriteVals(i);
                end
                tOverwrite(t,4) = toc;
            end

            % ============================================================
            %  5. KEY-EXISTENCE CHECK / KNOWN
            % ============================================================

            % --- struct (isfield) ---
            tic;
            for i = lookupOrder
                dummy = isfield(S, keys{i});
            end
            tHasKey(t,1) = toc;

            % --- containers.Map (isKey) ---
            tic;
            for i = lookupOrder
                dummy = M.isKey(keys{i});
            end
            tHasKey(t,2) = toc;

            % --- java.util.HashMap (containsKey) ---
            if testJava
                tic;
                for i = lookupOrder
                    dummy = J.containsKey(keys{i});
                end
                tHasKey(t,3) = toc;
            end

            % --- dictionary (isKey) ---
            if testDict
                tic;
                for i = lookupOrder
                    dummy = isKey(D, keys{i});
                end
                tHasKey(t,4) = toc;
            end

            % ============================================================
            %  6. KEY-EXISTENCE CHECK / UNKNOWN
            % ============================================================

            % --- struct (isfield) ---
            tic;
            for i = lookupOrder
                dummy = isfield(S, keys2{i});
            end
            tHasKey2(t,1) = toc;

            % --- containers.Map (isKey) ---
            tic;
            for i = lookupOrder
                dummy = M.isKey(keys2{i});
            end
            tHasKey2(t,2) = toc;

            % --- java.util.HashMap (containsKey) ---
            if testJava
                tic;
                for i = lookupOrder
                    dummy = J.containsKey(keys2{i});
                end
                tHasKey2(t,3) = toc;
            end

            % --- dictionary (isKey) ---
            if testDict
                tic;
                for i = lookupOrder
                    dummy = isKey(D, keys2{i});
                end
                tHasKey2(t,4) = toc;
            end

            % ============================================================
            %  7. ITERATION (visit every key-value pair)
            % ============================================================

            % --- struct ---
            tic;
            fn = fieldnames(S);
            acc = 0;
            for i = 1:numel(fn)
                acc = acc + S.(fn{i});
            end
            tIterate(t,1) = toc;

            % --- containers.Map ---
            tic;
            mk = M.keys;
            acc = 0;
            for i = 1:numel(mk)
                acc = acc + M(mk{i});
            end
            tIterate(t,2) = toc;

            % --- java.util.HashMap ---
            if testJava
                tic;
                it = J.entrySet.iterator;
                acc = 0;
                while it.hasNext
                    entry = it.next;
                    acc = acc + entry.getValue;
                end
                tIterate(t,3) = toc;
            end

            % --- dictionary ---
            if testDict
                tic;
                dk = D.keys();
                acc = 0;
                for i = 1:numel(dk)
                    acc = acc + D(dk(i));
                end
                tIterate(t,4) = toc;
            end

            % ============================================================
            %  8. DELETION
            % ============================================================
            deleteOrder = randperm(N);

            % --- struct (rmfield one-by-one) ---
            tic;
            for i = deleteOrder
                %S = rmfield(S, keys{i}); % this has a catastrophic performance!
                S.(keys{i})=[]; % this has implications for isKey() where isEmpty() nned to be added and memory freeing...
            end
            tDelete(t,1) = toc;

            % --- containers.Map (remove) ---
            tic;
            for i = deleteOrder
                M.remove(keys{i});
            end
            tDelete(t,2) = toc;

            % --- java.util.HashMap (remove) ---
            if testJava
                tic;
                for i = deleteOrder
                    J.remove(keys{i});
                end
                tDelete(t,3) = toc;
            end

            % --- dictionary (remove) ---
            if testDict
                tic; 
                for i = deleteOrder
                    D(keys{i})=[]; % D.remove(keys{i}); is catastrophically slow, similar to rmfield() for structs... but in contrast, D(keys{i})=[] actually removes the key, as isKey() returns then false 
                end
                tDelete(t,4) = toc;
            end
        end

        % Store median times
        row = [N, ...
               median(tInsert(:,1)),   median(tInsert(:,2)),   median(tInsert(:,3)),   median(tInsert(:,4)), ...
               median(tLookup(:,1)),   median(tLookup(:,2)),   median(tLookup(:,3)),   median(tLookup(:,4)), ...
               median(tLookup2(:,1)),  median(tLookup2(:,2)),  median(tLookup2(:,3)),  median(tLookup2(:,4)), ...
               median(tOverwrite(:,1)),median(tOverwrite(:,2)),median(tOverwrite(:,3)),median(tOverwrite(:,4)), ...
               median(tHasKey(:,1)),   median(tHasKey(:,2)),   median(tHasKey(:,3)),   median(tHasKey(:,4)), ...
               median(tHasKey2(:,1)),  median(tHasKey2(:,2)),  median(tHasKey2(:,3)),  median(tHasKey2(:,4)), ...
               median(tIterate(:,1)),  median(tIterate(:,2)),  median(tIterate(:,3)),  median(tIterate(:,4)), ...
               median(tDelete(:,1)),   median(tDelete(:,2)),   median(tDelete(:,3)),   median(tDelete(:,4))];
        results(si,:) = row;

        % Console summary for this size
        printBlock('Insert',    N, tInsert);
        printBlock('Lookup',    N, tLookup);
        printBlock('Lookup2',   N, tLookup2);
        printBlock('Overwrite', N, tOverwrite);
        printBlock('HasKey',    N, tHasKey);
        printBlock('HasKey2',   N, tHasKey2);
        printBlock('Iterate',   N, tIterate);
        printBlock('Delete',    N, tDelete);
    end

    % ==================================================================
    %  Final summary table
    % ==================================================================
    fprintf('\n\n==================== SUMMARY (median seconds) ====================\n');
    fprintf('%-6s | %-48s | %-48s | %-48s | %-48s | %-48s | %-48s | %-48s | %-48s\n', ...
        'N', 'Insert (S / M / J / D)', 'Lookup (S / M / J / D)', 'Lookup2 (S / M / J / D)', ...
        'Overwrite (S / M / J / D)', 'HasKey (S / M / J / D)', 'HasKey2 (S / M / J / D)', ...
        'Iterate (S / M / J / D)', 'Delete (S / M / J / D)');
    fprintf('%s\n', repmat('-', 1, 320));
    for si = 1:numel(sizes)
        r = results(si,:);
        fprintf('%-6d', r(1));
        for op = 0:7
            c = 2 + op*nImpl;
            fprintf(' | %10.5f / %10.5f / %10.5f / %10.5f', r(c), r(c+1), r(c+2), r(c+3));
        end
        fprintf('\n');
    end

    % ==================================================================
    %  Plot results
    % ==================================================================
    plotResults(sizes, results, nImpl);

    % ==================================================================
    %  Additional sweeps for N = 5000
    % ==================================================================
    benchmarkKeyLength(5000, nTrials, testJava, testDict);
    benchmarkPayloadSize(5000, nTrials, testJava, testDict);
end

% ======================================================================
%  KEY LENGTH SWEEP  (N fixed, vary key string length)
% ======================================================================
function benchmarkKeyLength(N, nTrials, testJava, testDict)
    keyLengths = [4, 8, 16, 32, 48, 63];
    % Note: MATLAB truncates struct field names to namelengthmax (63) chars,
    % so we cap at 63 to keep the comparison fair.
    % Operations: Insert, Lookup, Overwrite, HasKey
    % (skip Delete for struct – too slow at N=5000; skip Iterate – not key-sensitive)
    opNames = {'Insert','Lookup','Overwrite','HasKey'};
    nOps = numel(opNames);
    nImpl = 4;
    % results: rows = keyLengths, cols = [keyLen, S_ins M_ins J_ins D_ins, S_lkp ...]
    res = zeros(numel(keyLengths), 1 + nOps*nImpl);

    fprintf('\n\n############################################################\n');
    fprintf('  KEY LENGTH SWEEP  (N = %d, %d trials)\n', N, nTrials);
    fprintf('############################################################\n');

    for ki = 1:numel(keyLengths)
        kl = keyLengths(ki);
        fprintf('\n--- keyLen = %d ---\n', kl);
        keys = generateRandomKeys(N, kl);
        vals = randi(1e6, 1, N);
        overwriteVals = randi(1e6, 1, N);
        lookupOrder = randperm(N);

        tIns = zeros(nTrials, nImpl);
        tLkp = zeros(nTrials, nImpl);
        tOvr = zeros(nTrials, nImpl);
        tHas = zeros(nTrials, nImpl);

        for t = 1:nTrials
            % --- Insert ---
            S = struct(); tic;
            for i = 1:N, S.(keys{i}) = vals(i); end
            tIns(t,1) = toc;

            tic; M = containers.Map(keys, vals); tIns(t,2) = toc;

            if testJava
                J = java.util.HashMap(N); tic;
                for i = 1:N, J.put(keys{i}, vals(i)); end
                tIns(t,3) = toc;
            end

            if testDict
                tic; D = dictionary(string(keys), vals); tIns(t,4) = toc;
            end

            % --- Lookup ---
            dummy = 0; %#ok<NASGU>
            tic; for i = lookupOrder, dummy = S.(keys{i}); end; tLkp(t,1) = toc;
            tic; for i = lookupOrder, dummy = M(keys{i});   end; tLkp(t,2) = toc;
            if testJava
                tic; for i = lookupOrder, dummy = J.get(keys{i});end; tLkp(t,3) = toc;
            end
            if testDict
                tic; for i = lookupOrder, dummy = D(keys{i}); end; tLkp(t,4) = toc;
            end

            % --- Overwrite ---
            tic; for i = 1:N, S.(keys{i}) = overwriteVals(i); end; tOvr(t,1) = toc;
            tic; for i = 1:N, M(keys{i})  = overwriteVals(i); end; tOvr(t,2) = toc;
            if testJava
                tic; for i = 1:N, J.put(keys{i}, overwriteVals(i)); end; tOvr(t,3) = toc;
            end
            if testDict
                tic; for i = 1:N, D(keys{i}) = overwriteVals(i); end; tOvr(t,4) = toc;
            end

            % --- HasKey ---
            tic; for i = lookupOrder, dummy = isfield(S, keys{i});      end; tHas(t,1) = toc;
            tic; for i = lookupOrder, dummy = M.isKey(keys{i});         end; tHas(t,2) = toc;
            if testJava
                tic; for i = lookupOrder, dummy = J.containsKey(keys{i});   end; tHas(t,3) = toc;
            end
            if testDict
                tic; for i = lookupOrder, dummy = isKey(D, keys{i}); end; tHas(t,4) = toc;
            end
        end

        res(ki,:) = [kl, ...
            median(tIns(:,1)), median(tIns(:,2)), median(tIns(:,3)), median(tIns(:,4)), ...
            median(tLkp(:,1)), median(tLkp(:,2)), median(tLkp(:,3)), median(tLkp(:,4)), ...
            median(tOvr(:,1)), median(tOvr(:,2)), median(tOvr(:,3)), median(tOvr(:,4)), ...
            median(tHas(:,1)), median(tHas(:,2)), median(tHas(:,3)), median(tHas(:,4))];

        printBlock('Insert',   kl, tIns);
        printBlock('Lookup',   kl, tLkp);
        printBlock('Overwrite',kl, tOvr);
        printBlock('HasKey',   kl, tHas);
    end

    % --- table ---
    fprintf('\n  KEY LENGTH SWEEP SUMMARY (N=%d, median seconds)\n', N);
    fprintf('  %-8s', 'KeyLen');
    for o = 1:nOps
        fprintf(' | %-46s', [opNames{o} ' (S / M / J / D)']);
    end
    fprintf('\n  %s\n', repmat('-', 1, 8 + nOps*49));
    for ki = 1:numel(keyLengths)
        r = res(ki,:);
        fprintf('  %-8d', r(1));
        for o = 1:nOps
            c = 2 + (o-1)*nImpl;
            fprintf(' | %10.5f / %10.5f / %10.5f / %10.5f', r(c), r(c+1), r(c+2), r(c+3));
        end
        fprintf('\n');
    end

    % --- plot ---
    if mr.aux.isOctave()
        add='_octave';
    else
        add='_matlab';
    end
    plotSweep(keyLengths, res, opNames, nImpl, ...
        'Key Length (chars)', ...
        sprintf('Key Length Sweep (N=%d)', N), ...
        ['benchmarkKeyLength' add '.png']);
end

% ======================================================================
%  PAYLOAD SIZE SWEEP  (N fixed, vary value vector length)
% ======================================================================
function benchmarkPayloadSize(N, nTrials, testJava, testDict)
    payloadSizes = [1, 10, 100, 1000, 10000];
    keyLen = 12;
    opNames = {'Insert','Lookup','Overwrite','Iterate'};
    nOps = numel(opNames);
    nImpl = 4;
    res = zeros(numel(payloadSizes), 1 + nOps*nImpl);

    fprintf('\n\n############################################################\n');
    fprintf('  PAYLOAD SIZE SWEEP  (N = %d, keyLen = %d, %d trials)\n', N, keyLen, nTrials);
    fprintf('############################################################\n');

    keys = generateRandomKeys(N, keyLen);

    for pi = 1:numel(payloadSizes)
        pLen = payloadSizes(pi);
        fprintf('\n--- payload length = %d ---\n', pLen);

        % Generate payload values: each value is a 1xpLen double vector
        vals = cell(1, N);
        overwriteVals = cell(1, N);
        for i = 1:N
            vals{i} = randi(1e6, 1, pLen);
            overwriteVals{i} = randi(1e6, 1, pLen);
        end
        lookupOrder = randperm(N);

        tIns = zeros(nTrials, nImpl);
        tLkp = zeros(nTrials, nImpl);
        tOvr = zeros(nTrials, nImpl);
        tItr = zeros(nTrials, nImpl);

        for t = 1:nTrials
            % ---- Insert ----
            S = struct(); tic;
            for i = 1:N, S.(keys{i}) = vals{i}; end
            tIns(t,1) = toc;

            tic;
            M = containers.Map();
            for i = 1:N, M(keys{i}) = vals{i}; end
            tIns(t,2) = toc;

            if testJava
                J = java.util.HashMap(N); tic;
                for i = 1:N, J.put(keys{i}, vals{i}); end
                tIns(t,3) = toc;
            end

            if testDict
                tic;
                D = configureDictionary("string","cell");
                for i = 1:N, D(keys{i}) = {vals{i}}; end
                tIns(t,4) = toc;
            end

            % ---- Lookup ----
            dummy = 0; %#ok<NASGU>
            tic; for i = lookupOrder, dummy = S.(keys{i}); end; tLkp(t,1) = toc;
            tic; for i = lookupOrder, dummy = M(keys{i});  end; tLkp(t,2) = toc;
            if testJava
                tic; for i = lookupOrder, dummy = J.get(keys{i}); end; tLkp(t,3) = toc;
            end
            if testDict
                tic; for i = lookupOrder, dummy = D(keys{i}); end; tLkp(t,4) = toc;
            end

            % ---- Overwrite ----
            tic; for i = 1:N, S.(keys{i}) = overwriteVals{i}; end; tOvr(t,1) = toc;
            tic; for i = 1:N, M(keys{i})  = overwriteVals{i}; end; tOvr(t,2) = toc;
            if testJava
                tic; for i = 1:N, J.put(keys{i}, overwriteVals{i}); end; tOvr(t,3) = toc;
            end
            if testDict
                tic; for i = 1:N, D(keys{i}) = {overwriteVals{i}}; end; tOvr(t,4) = toc;
            end

            % ---- Iterate ----
            tic;
            fn = fieldnames(S); acc = 0;
            for i = 1:numel(fn), acc = acc + sum(S.(fn{i})); end
            tItr(t,1) = toc;

            tic;
            mk = M.keys; acc = 0;
            for i = 1:numel(mk), acc = acc + sum(M(mk{i})); end
            tItr(t,2) = toc;

            if testJava
                tic;
                it = J.entrySet.iterator; acc = 0;
                while it.hasNext
                    entry = it.next;
                    v = double(entry.getValue);
                    acc = acc + sum(v);
                end
                tItr(t,3) = toc;
            end

            if testDict
                tic;
                dk = D.keys(); acc = 0;
                for i = 1:numel(dk)
                    v = D(dk(i));
                    if iscell(v), v = v{1}; end
                    acc = acc + sum(v);
                end
                tItr(t,4) = toc;
            end
        end

        res(pi,:) = [pLen, ...
            median(tIns(:,1)), median(tIns(:,2)), median(tIns(:,3)), median(tIns(:,4)), ...
            median(tLkp(:,1)), median(tLkp(:,2)), median(tLkp(:,3)), median(tLkp(:,4)), ...
            median(tOvr(:,1)), median(tOvr(:,2)), median(tOvr(:,3)), median(tOvr(:,4)), ...
            median(tItr(:,1)), median(tItr(:,2)), median(tItr(:,3)), median(tItr(:,4))];

        printBlock('Insert',   pLen, tIns);
        printBlock('Lookup',   pLen, tLkp);
        printBlock('Overwrite',pLen, tOvr);
        printBlock('Iterate',  pLen, tItr);
    end

    % --- table ---
    fprintf('\n  PAYLOAD SIZE SWEEP SUMMARY (N=%d, median seconds)\n', N);
    fprintf('  %-10s', 'Payload');
    for o = 1:nOps
        fprintf(' | %-46s', [opNames{o} ' (S / M / J / D)']);
    end
    fprintf('\n  %s\n', repmat('-', 1, 10 + nOps*49));
    for pi = 1:numel(payloadSizes)
        r = res(pi,:);
        fprintf('  %-10d', r(1));
        for o = 1:nOps
            c = 2 + (o-1)*nImpl;
            fprintf(' | %10.5f / %10.5f / %10.5f / %10.5f', r(c), r(c+1), r(c+2), r(c+3));
        end
        fprintf('\n');
    end

    % --- plot ---
    if mr.aux.isOctave()
        add='_octave';
    else
        add='_matlab';
    end
    plotSweep(payloadSizes, res, opNames, nImpl, ...
        'Payload Size (doubles)', ...
        sprintf('Payload Size Sweep (N=%d)', N), ...
        ['benchmarkPayloadSize' add '.png']);
end

% ======================================================================
%  Helper: detect dictionary availability (R2022b+)
% ======================================================================
function tf = hasDictionary()
    tf = false;
    try
        d = dictionary(string("a"), 1); %#ok<NASGU>
        tf = true;
    catch
    end
end

% ======================================================================
%  Helper: plot a parameter-sweep figure (reused by both sweeps)
% ======================================================================
function plotSweep(xVals, res, opNames, nImpl, xLabel, figTitle, fileName)
    nOps = numel(opNames);
    nCols = min(nOps, 3);
    nRows = ceil(nOps / nCols);
    legLabels = {'struct','containers.Map','Java HashMap','dictionary'};
    colors = [0.2 0.6 0.9; 0.9 0.4 0.1; 0.3 0.8 0.3; 0.8 0.2 0.8];

    figure('Name', figTitle, 'Position', [100 100 400*nCols 350*nRows]);
    for p = 1:nOps
        subplot(nRows, nCols, p);
        c = 2 + (p-1)*nImpl;
        data = res(:, c:c+nImpl-1);
        cats = arrayfun(@num2str, xVals, 'UniformOutput', false);
        if mr.aux.isOctave()
          b = bar(data);
          set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
        else
          b = bar(cats, data);
          for bi = 1:numel(b)
              b(bi).FaceColor = colors(bi,:);
          end
        end
        ylabel('Time (s)');
        xlabel(xLabel);
        title(opNames{p});
        legend(legLabels, 'Location','northwest');
        grid on;
    end
    if ~mr.aux.isOctave()
      sgtitle(figTitle);
    end
    drawnow;
    try
        saveas(gcf, fileName);
        fprintf('\nFigure saved to %s\n', fileName);
    catch
        fprintf('\n(Could not save figure – possibly headless environment)\n');
    end
end

% ======================================================================
%  Helper: generate N unique random strings that are valid identifiers
% ======================================================================
function keys = generateRandomKeys(N, keyLen)
    chars = ['a':'z', 'A':'Z'];
    keys  = cell(1, N);
    used  = containers.Map('KeyType','char','ValueType','logical');
    for i = 1:N
        while true
            k = chars(randi(numel(chars), 1, keyLen));
            if ~used.isKey(k)
                used(k) = true;
                keys{i} = k;
                break;
            end
        end
    end
end

% ======================================================================
%  Helper: print a timing block
% ======================================================================
function printBlock(label, N, tMatrix)
    med = median(tMatrix, 1);
    nCol = size(tMatrix, 2);
    if nCol >= 4
        fprintf('  %-10s (N=%5d):  struct %8.5f s  |  Map %8.5f s  |  Java %8.5f s  |  Dict %8.5f s\n', ...
            label, N, med(1), med(2), med(3), med(4));
    else
        fprintf('  %-10s (N=%5d):  struct %8.5f s  |  Map %8.5f s  |  Java %8.5f s\n', ...
            label, N, med(1), med(2), med(3));
    end
end

% ======================================================================
%  Helper: plot comparative bar charts
% ======================================================================
function plotResults(sizes, results, nImpl)
    opNames = {'Insert','Lookup','Lookup2','Overwrite','HasKey','HasKey2','Iterate','Delete'};
    nOps = numel(opNames);
    % Build column indices dynamically: col 1 is N, then nImpl cols per op
    cols = cell(1, nOps);
    for p = 1:nOps
        c = 2 + (p-1)*nImpl;
        cols{p} = c : c+nImpl-1;
    end
    legLabels = {'struct','containers.Map','Java HashMap','dictionary'};
    colors = [0.2 0.6 0.9; 0.9 0.4 0.1; 0.3 0.8 0.3; 0.8 0.2 0.8];

    figure('Name','Map Performance Comparison','Position',[100 100 1800 800]);

    for p = 1:nOps
        subplot(2, 4, p);
        data = results(:, cols{p});
        cats = arrayfun(@num2str, sizes, 'UniformOutput', false);
        if mr.aux.isOctave()
          b = bar(data);
          set(gca,'XTick',1:numel(cats),'XTickLabel',cats);
        else
          b = bar(cats, data);
          for bi = 1:numel(b)
              b(bi).FaceColor = colors(bi,:);
          end
        end
        ylabel('Time (s)');
        xlabel('N (number of keys)');
        title(opNames{p});
        legend(legLabels, 'Location','northwest');
        grid on;
    end

    if ~mr.aux.isOctave()
      sgtitle('Struct vs containers.Map vs Java HashMap vs dictionary – Benchmark');
    end
    drawnow;

    if mr.aux.isOctave()
        add='_octave';
    else
        add='_matlab';
    end

    % Also save to file for headless runs
    try
        saveas(gcf, ['benchmarkMapPerformance' add '.png']);
        fprintf('\nFigure saved to benchmarkMapPerformance%s.png\n',add);
    catch
        fprintf('\n(Could not save figure – possibly headless environment)\n');
    end
end
