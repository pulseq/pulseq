classdef testRotationExtension < matlab.unittest.TestCase
    % testRotationExtension
    %   This class implements tests for a radial sequence, its rotations,
    %   write/read, plotting, and recreation. It is a 1–1 translation of your
    %   Python test suite.
    
    properties (Constant)
        % Expected output directory (relative to this file)
        expectedOutputPath = fullfile(fileparts(mfilename('fullpath')), 'expected_output');
    end
    
    methods (Test)
        
        function test_vs_rotate(testCase)
            % Test that explicit gradient rotation (via pp.rotate) yields the same results
            % as the version using rotation extensions.
            seq  = seq_make_radial();
            seq2 = seq_make_radial_norotext();
            
            % Test waveforms_and_times (assumed to return cell arrays of cells)
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
                       
             % Test approximate equality of k-space calculation.
            kspace1 = seq.calculateKspacePP();
            kspace2 = seq2.calculateKspacePP();
            verifyApproxEqual(testCase, kspace2, kspace1, 1e-1, []);
        end
        
        function test_sequence_save_expected(testCase)
            % When the SAVE_EXPECTED environment variable is set,
            % write the sequence file to the expected output directory.
            if isempty(getenv('SAVE_EXPECTED'))
                testCase.assumeFail('SAVE_EXPECTED not set; skipping saving expected sequence files.');
            end
            seq = seq_make_radial();
            filename = fullfile(testRotationExtension.expectedOutputPath, 'seq_make_radial.seq');
            seq.write(filename);
            % (Additional verification could be added here.)
        end
        
        function test_plot(testCase)
            % Test that the sequence can be plotted without errors.
            seq = seq_make_radial();
            try
                seq.plot();
                seq.plot('showBlocks', true);
                seq.plot('timeRange', [0, 1e-3]);
                seq.plot('timeDisp', 'ms');
                close all;
            catch ME
                testCase.verifyFail(['Plotting failed: ', ME.message]);
            end
        end
        
        function test_writeread(testCase)
            % Test that writing a sequence to file and reading it back
            % produces an equivalent sequence.
            tempFile = fullfile(tempdir, 'seq_make_radial.seq');
            seq = seq_make_radial();
            
            % Write sequence to file
            seq.write(tempFile);
            
            % Read written sequence back in
            seq2 = mr.Sequence();  % New sequence with same system settings.
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
        
        function test_recreate(testCase)
            % Test that recreating a sequence by extracting its blocks and re‑adding them
            % to a new sequence yields an equivalent sequence.
            seq = seq_make_radial();
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
end

function R = rotation_matrix(angle)
    % ROTATION_MATRIX Create a 3x3 rotation matrix for a given angle (degrees)
    theta = deg2rad(angle);
    R0 = [cos(theta), -sin(theta), 0];
    R1 = [sin(theta),  cos(theta), 0];
    R2 = [0, 0, 1];
    R = [R0; R1; R2];
end

function seq = seq_make_radial()
    % SEQ_MAKE_RADIAL creates a basic radial sequence.
    seq = mr.Sequence();
    rf = mr.makeBlockPulse(pi/2, 'duration', 1e-3);
    gread = mr.makeTrapezoid('x', 'area', 1000);
    theta = [0, 30, 45, 60, 90];
    rot = cell(1, numel(theta));
    for n = 1:numel(theta)
        rot{n} = rotation_matrix(theta(n));
    end
    for n = 1:numel(theta)
        seq.addBlock(rf);
        seq.addBlock(gread, mr.makeRotation(rot{n}));
    end
    seq.addBlock(rf);
    seq.addBlock(gread, mr.makeRotation(rot{1}));
end

function seq = seq_make_radial_norotext()
    % SEQ_MAKE_RADIAL_NOROTEXT creates a radial sequence using explicit rotation
    % via the pp.rotate function.
    seq = mr.Sequence();
    rf = mr.makeBlockPulse(pi/2, 'duration', 1e-3);
    gread = mr.makeTrapezoid('x', 'area', 1000);
    
    theta = [0, 30, 45, 60, 90];
    theta_rad = deg2rad(theta);
    for n = 1:numel(theta_rad)
        seq.addBlock(rf);
        
        args = mr.rotate('z', theta_rad(n), gread);
        seq.addBlock(args{:});
    end
    seq.addBlock(rf);
    
    args = mr.rotate('z', theta_rad(1), gread);
    seq.addBlock(args{:});
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
