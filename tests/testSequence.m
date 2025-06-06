classdef testSequence < matlab.unittest.TestCase
    % TestSequence
    %   This class implements tests for sequence generation, writing/reading,
    %   plotting, and recreating the sequence, analogous to the original Python tests.
    %
    %   It uses a parameterized property "seqFunc" that loops over a zoo of sequence‐
    %   creation functions. Example sequences (from Python scripts) are also added.
    
    properties (TestParameter)
        % TestParameter "seqFunc" will take on each sequence‐creation function
        seqFunc = testSequence.getSequenceZoo();
    end
    
    methods (Test)
        function test_sequence(testCase, seqFunc)
            % Base test that runs the sequence function and verifies the result is nonempty.
            seq = seqFunc();
            testCase.verifyNotEmpty(seq, 'Sequence is empty.');
            % Store the sequence for later tests (if needed)
            assignin('base', 'TestSequence_seq', seq);
        end
        
        function test_save_expected(testCase, seqFunc)
            % This test rewrites the expected .seq output file when SAVE_EXPECTED is set.
            if isempty(getenv('SAVE_EXPECTED'))
                testCase.assumeFail('SAVE_EXPECTED not set; skipping saving expected sequence files.');
            end
            seq = seqFunc();
            seq_name = func2str(seqFunc);
            % Determine expected output file path relative to this file.
            testFileDir = fileparts(mfilename('fullpath'));
            expectedDir = fullfile(testFileDir, 'expected_output');
            if ~exist(expectedDir, 'dir')
                mkdir(expectedDir);
            end
            expected_file = fullfile(expectedDir, [seq_name, '.seq']);
            seq.write(expected_file);
        end
        
        function test_plot(testCase, seqFunc)
            % Test the sequence.plot() method for selected sequences.
            seq_name = func2str(seqFunc);
            % Only run for these sequences.
            if any(strcmp(seq_name, {'seq1', 'seq2', 'seq3', 'seq4'}))
                seq = seqFunc();
                try
                    % These plotting options should not error.
                    seq.plot();
                    seq.plot('showBlocks', true);
                    seq.plot('timeRange', [0, 1e-3]);
                    seq.plot('timeDisp', 'ms');
                    close all;
                catch ME
                    testCase.verifyFail(['Plotting failed: ', ME.message]);
                end
            else
                testCase.assumeFail('Plot test not applicable for this sequence.');
            end
            close all
        end
        
        function test_writeread(testCase, seqFunc)
            % Test whether the sequence remains approximately the same after writing and reading.
            seq_name = func2str(seqFunc);
            tempFile = fullfile(tempdir, [seq_name, '.seq']);
            seq = seqFunc();
            
            % Write sequence to file
            seq.write(tempFile);
            
            % Read written sequence back in
            seq2 = mr.Sequence();
            seq2.read(tempFile);
            
            % Compare sequence block events (using field names as keys).
            testCase.verifyEqual(length(seq.blockEvents), length(seq2.blockEvents), ...
                'Sequence block IDs are not identical');
            
            for i = 1:length(seq.blockEvents)
                block_orig = seq.getBlock(i);
                block_compare = seq2.getBlock(i);
                
                % If block contains rf.use, set it to 'undefined'
                if isfield(block_orig, 'rf') && isfield(block_orig.rf, 'use')
                    block_orig.rf.use = 'undefined';
                end
                
                verifyApproxEqual(testCase, block_compare, block_orig, 1e-5, 1e-5);
            end
            
            % Compare gradient waveforms.
            grads1 = seq.waveforms_and_times();
            grads2 = seq2.waveforms_and_times();
            channels = {'x', 'y', 'z'};
            for ch = 1:length(channels)
                if isempty(grads1{ch}) && isempty(grads2{ch})
                    continue;
                end
                verifyApproxEqual(testCase, grads1{ch}(1, :), grads2{ch}(1, :), 1e-5, 1e-5, ...
                    sprintf('Time axis of gradient waveform for channel %s does not match', channels{ch}));
                verifyApproxEqual(testCase, grads1{ch}(2, :), grads2{ch}(2, :), 1e2, 1e-3, ...
                    sprintf('Gradient values of gradient waveform for channel %s do not match', channels{ch}));
            end
            
            % Restore RF use for k-space calculation.
            for i = 1:length(seq.blockEvents)
                block_orig = seq.getBlock(i);
                if isfield(block_orig, 'rf') && isfield(block_orig.rf, 'use')
                    block_compare = seq2.getBlock(i);
                    block_compare.rf.use = block_orig.rf.use;
                end
            end
            
            % Test approximate equality of k-space calculation.
            kspace1 = seq.calculateKspacePP();
            kspace2 = seq2.calculateKspacePP();
            verifyApproxEqual(testCase, kspace2, kspace1, 1e-1, []);
            
            % Test whether labels are the same.
            labels_seq = seq.evalLabels('evolution', 'blocks');
            labels_seq2 = seq2.evalLabels('evolution', 'blocks');
            testCase.verifyEqual(fieldnames(labels_seq), fieldnames(labels_seq2), 'Sequences do not contain the same set of labels');
            label_keys = fieldnames(labels_seq);
            for i = 1:length(label_keys)
                key = label_keys{i};
                testCase.verifyEqual(labels_seq.(key), labels_seq2.(key), ...
                    sprintf('Label %s does not match', key));
            end
            
            if exist(tempFile, 'file')
                delete(tempFile);
            end
        end
        
        function test_recreate(testCase, seqFunc)
            % Test whether the sequence is approximately the same after recreating
            % it by getting all blocks and inserting them into a new sequence.
            seq = seqFunc();
            seq2 = mr.Sequence();
            for i = 1:length(seq.blockEvents)
                seq2.addBlock(seq.getBlock(i));
            end
            for i = 1:length(seq.blockEvents)
                verifyApproxEqual(testCase, seq2.getBlock(i), seq.getBlock(i), 1e-9, 1e-9, ...
                    sprintf('Block %s does not match', i));
            end
            
            % Compare gradient waveforms.
            grads1 = seq.waveforms_and_times();
            grads2 = seq2.waveforms_and_times();
            channels = {'x','y','z'};
            for ch = 1:length(channels)
                if isempty(grads1{ch}) && isempty(grads2{ch})
                    continue;
                end
                verifyApproxEqual(testCase, grads1{ch}(1, :), grads2{ch}(1, :), 1e-9, 1e-9, ...
                    sprintf('Time axis of gradient waveform for channel %s does not match', channels{ch}));
                verifyApproxEqual(testCase, grads1{ch}(2, :), grads2{ch}(2, :), 1e-9, 1e-9, ...
                    sprintf('Gradient values of gradient waveform for channel %s do not match', channels{ch}));
            end
            
            verifyApproxEqual(testCase, seq2.calculateKspacePP(), seq.calculateKspacePP(), 1e-6, []);
            
            labels_seq = seq.evalLabels('evolution', 'blocks');
            labels_seq2 = seq2.evalLabels('evolution', 'blocks');
            testCase.verifyEqual(fieldnames(labels_seq), fieldnames(labels_seq2), 'Sequences do not contain the same set of labels');
            label_keys = fieldnames(labels_seq);
            for i = 1:length(label_keys)
                key = label_keys{i};
                testCase.verifyEqual(labels_seq.(key), labels_seq2.(key), ...
                    sprintf('Label %s does not match', key));
            end
        end
    end
    
    methods (Static)
        function zoo = getSequenceZoo()
            % Returns a cell array of function handles that create sequences.
            zoo = {@seq_make_gauss_pulses, ...
                   @seq_make_sinc_pulses, ...
                   @seq_make_block_pulses, ...
                   @seq1, ...
                   @seq2, ...
                   @seq3, ...
                   @seq4 ...
                  };
        end
        
    end
end

%% --- Auxiliary Functions ---
function seq = seq_make_gauss_pulses()
    seq = mr.Sequence();
    seq.addBlock(mr.makeGaussPulse(1, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(1, 'duration', 1e-3, 'delay', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(pi/2, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(pi/2, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(pi/2, 'duration', 2e-3, 'phaseOffset', pi/2));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(pi/2, 'duration', 1e-3, 'phaseOffset', pi/2, 'freqOffset', 1e3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(pi/2, 'duration', 1e-3, 'timeBwProduct', 1));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeGaussPulse(pi/2, 'duration', 1e-3, 'apodization', 0.1));
end

function seq = seq_make_sinc_pulses()
    seq = mr.Sequence();
    seq.addBlock(mr.makeSincPulse(1, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(1, 'duration', 1e-3, 'delay', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(pi/2, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(pi/2, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(pi/2, 'duration', 2e-3, 'phaseOffset', pi/2));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(pi/2, 'duration', 1e-3, 'phaseOffset', pi/2, 'freqOffset', 1e3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(pi/2, 'duration', 1e-3, 'timeBwProduct', 1));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeSincPulse(pi/2, 'duration', 1e-3, 'apodization', 0.1));
end

function seq = seq_make_block_pulses()
    seq = mr.Sequence();
    seq.addBlock(mr.makeBlockPulse(1, 'duration', 4e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeBlockPulse(1, 'duration', 4e-3, 'delay', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeBlockPulse(pi/2, 'duration', 4e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeBlockPulse(pi/2, 'duration', 1e-3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeBlockPulse(pi/2, 'duration', 2e-3, 'phaseOffset', pi/2));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeBlockPulse(pi/2, 'duration', 1e-3, 'phaseOffset', pi/2, 'freqOffset', 1e3));
    seq.addBlock(mr.makeDelay(1));
    seq.addBlock(mr.makeBlockPulse(pi/2, 'duration', 1e-3, 'timeBwProduct', 1));
end

function seq = seq1()
    % Basic sequence with gradients in all channels.
    seq = mr.Sequence();
    seq.addBlock(mr.makeBlockPulse(pi/4, 'duration', 1e-3));
    seq.addBlock(mr.makeTrapezoid('x', 'area', 1000));
    seq.addBlock(mr.makeTrapezoid('y', 'area', -500.00001));
    seq.addBlock(mr.makeTrapezoid('z', 'area', 100));
    seq.addBlock(mr.makeTrapezoid('x', 'area', -1000), mr.makeTrapezoid('y', 'area', 500));
    seq.addBlock(mr.makeTrapezoid('y', 'area', -500), mr.makeTrapezoid('z', 'area', 1000));
    seq.addBlock(mr.makeTrapezoid('x', 'area', -1000), mr.makeTrapezoid('z', 'area', 1000.00001));
end

function seq = seq2()
    % Basic spin-echo sequence structure.
    seq = mr.Sequence();
    seq.addBlock(mr.makeBlockPulse(pi/2, 'duration', 1e-3));
    seq.addBlock(mr.makeTrapezoid('x', 'area', 1000));
    seq.addBlock(mr.makeTrapezoid('x', 'area', -1000));
    seq.addBlock(mr.makeBlockPulse(pi, 'duration', 1e-3));
    seq.addBlock(mr.makeTrapezoid('x', 'area', -500));
    seq.addBlock(mr.makeTrapezoid('x', 'area', 1000, 'duration', 10e-3), ...
                 mr.makeAdc(100, 'duration', 10e-3));
end

function seq = seq3()
    % Basic GRE sequence with INC labels.
    seq = mr.Sequence();
    for i = 0:9
        seq.addBlock(mr.makeBlockPulse(pi/8, 'duration', 1e-3));
        seq.addBlock(mr.makeTrapezoid('x', 'area', 1000));
        seq.addBlock(mr.makeTrapezoid('y', 'area', -500 + i * 100));
        seq.addBlock(mr.makeTrapezoid('x', 'area', -500));
        seq.addBlock(mr.makeTrapezoid('x', 'area', 1000, 'duration', 10e-3), ...
                     mr.makeAdc(100, 'duration', 10e-3), ...
                     mr.makeLabel('INC', 'LIN', 1));
    end
end

function seq = seq4()
    % Basic GRE sequence with SET labels.
    seq = mr.Sequence();
    for i = 0:9
        seq.addBlock(mr.makeBlockPulse(pi/8, 'duration', 1e-3));
        seq.addBlock(mr.makeTrapezoid('x', 'area', 1000));
        seq.addBlock(mr.makeTrapezoid('y', 'area', -500 + i * 100));
        seq.addBlock(mr.makeTrapezoid('x', 'area', -500));
        seq.addBlock(mr.makeTrapezoid('x', 'area', 1000, 'duration', 10e-3), ...
                     mr.makeAdc(100, 'duration', 10e-3), ...
                     mr.makeLabel('SET', 'LIN', i));
    end
end

%% --- Helper Function for Approximate Comparison ---
function verifyApproxEqual(testCase, actual, expected, atol, rtol, varargin)
    % Recursively verify approximate equality between actual and expected.
    if isnumeric(actual) && isnumeric(expected)
        if isempty(rtol)
            testCase.verifyEqual(actual, expected, 'AbsTol', atol, varargin{:});
        else
            testCase.verifyEqual(actual, expected, 'AbsTol', atol, 'RelTol', rtol, varargin{:});
        end
    elseif isstruct(actual) && isstruct(expected)
        f1 = sort(fieldnames(actual));
        f2 = sort(fieldnames(expected));
        testCase.verifyEqual(f1, f2, 'Structures have different fields.');
        for i = 1:length(f1)
            fld = f1{i};
            verifyApproxEqual(testCase, actual.(fld), expected.(fld), atol, rtol, varargin{:});
        end
    elseif iscell(actual) && iscell(expected)
        testCase.verifyEqual(numel(actual), numel(expected), 'Cell arrays differ in length.');
        for i = 1:numel(actual)
            verifyApproxEqual(testCase, actual{i}, expected{i}, atol, rtol, varargin{:});
        end
    else
        if isempty(rtol)
            testCase.verifyEqual(actual, expected, 'AbsTol', atol, varargin{:});
        else
            testCase.verifyEqual(actual, expected, 'AbsTol', atol, 'RelTol', rtol, varargin{:});
        end
    end
end
