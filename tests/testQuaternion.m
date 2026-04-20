classdef testQuaternion < matlab.unittest.TestCase
    % Tests for mr.aux.quat.* functions:
    %   normalize, conjugate, multiply, toRotMat, fromRotMat, rotate

    methods (Test)

        %% normalize
        function test_normalize_unit(testCase)
            q = mr.aux.quat.normalize([1 0 0 0]);
            testCase.verifyEqual(norm(q), 1, 'AbsTol', 1e-12);
        end

        function test_normalize_arbitrary(testCase)
            q = mr.aux.quat.normalize([2 0 0 0]);
            testCase.verifyEqual(q, [1 0 0 0], 'AbsTol', 1e-12);
        end

        function test_normalize_general(testCase)
            q = mr.aux.quat.normalize([1 1 1 1]);
            testCase.verifyEqual(norm(q), 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(q, [0.5 0.5 0.5 0.5], 'AbsTol', 1e-12);
        end

        %% conjugate
        function test_conjugate(testCase)
            q = [1 2 3 4];
            qc = mr.aux.quat.conjugate(q);
            testCase.verifyEqual(qc, [1 -2 -3 -4]);
        end

        function test_conjugate_identity(testCase)
            q = [1 0 0 0];
            qc = mr.aux.quat.conjugate(q);
            testCase.verifyEqual(qc, [1 0 0 0]);
        end

        %% multiply
        function test_multiply_identity(testCase)
            q = mr.aux.quat.normalize([1 2 3 4]);
            id = [1 0 0 0];
            result = mr.aux.quat.multiply(q, id);
            testCase.verifyEqual(result, q, 'AbsTol', 1e-12);
        end

        function test_multiply_conjugate_gives_identity(testCase)
            q = mr.aux.quat.normalize([1 2 3 4]);
            qc = mr.aux.quat.conjugate(q);
            result = mr.aux.quat.multiply(q, qc);
            testCase.verifyEqual(result, [1 0 0 0], 'AbsTol', 1e-10);
        end

        function test_multiply_associativity(testCase)
            q1 = mr.aux.quat.normalize([1 1 0 0]);
            q2 = mr.aux.quat.normalize([1 0 1 0]);
            q3 = mr.aux.quat.normalize([1 0 0 1]);
            r1 = mr.aux.quat.multiply(mr.aux.quat.multiply(q1, q2), q3);
            r2 = mr.aux.quat.multiply(q1, mr.aux.quat.multiply(q2, q3));
            testCase.verifyEqual(r1, r2, 'AbsTol', 1e-10);
        end

        %% toRotMat
        function test_toRotMat_identity(testCase)
            q = [1 0 0 0];
            R = mr.aux.quat.toRotMat(q);
            testCase.verifyEqual(R, eye(3), 'AbsTol', 1e-12);
        end

        function test_toRotMat_90deg_z(testCase)
            % 90 degrees about z: q = [cos(pi/4), 0, 0, sin(pi/4)]
            q = [cos(pi/4), 0, 0, sin(pi/4)];
            R = mr.aux.quat.toRotMat(q);
            % Should rotate x -> y
            v_rotated = R * [1; 0; 0];
            testCase.verifyEqual(v_rotated, [0; 1; 0], 'AbsTol', 1e-10);
        end

        function test_toRotMat_180deg_x(testCase)
            % 180 degrees about x: q = [0, 1, 0, 0]
            q = [0, 1, 0, 0];
            R = mr.aux.quat.toRotMat(q);
            % Should flip y and z
            v = R * [0; 1; 0];
            testCase.verifyEqual(v, [0; -1; 0], 'AbsTol', 1e-10);
        end

        function test_toRotMat_orthogonal(testCase)
            q = mr.aux.quat.normalize([1 1 1 1]);
            R = mr.aux.quat.toRotMat(q);
            testCase.verifyEqual(R * R', eye(3), 'AbsTol', 1e-10);
            testCase.verifyEqual(det(R), 1, 'AbsTol', 1e-10);
        end

        %% fromRotMat
        function test_fromRotMat_identity(testCase)
            q = mr.aux.quat.fromRotMat(eye(3));
            testCase.verifyEqual(abs(q), [1 0 0 0], 'AbsTol', 1e-10);
        end

        function test_fromRotMat_roundtrip(testCase)
            q_orig = mr.aux.quat.normalize([1 2 3 4]);
            R = mr.aux.quat.toRotMat(q_orig);
            q_recov = mr.aux.quat.fromRotMat(R);
            % Quaternions q and -q represent same rotation
            if dot(q_orig, q_recov) < 0
                q_recov = -q_recov;
            end
            testCase.verifyEqual(q_recov, q_orig, 'AbsTol', 1e-10);
        end

        function test_fromRotMat_90deg_x(testCase)
            Rx = [1 0 0; 0 0 -1; 0 1 0];
            q = mr.aux.quat.fromRotMat(Rx);
            R_back = mr.aux.quat.toRotMat(q);
            testCase.verifyEqual(R_back, Rx, 'AbsTol', 1e-10);
        end

        %% rotate
        function test_rotate_identity(testCase)
            q = [1 0 0 0];
            v = [1 2 3];
            v_rot = mr.aux.quat.rotate(q, v);
            testCase.verifyEqual(v_rot, v, 'AbsTol', 1e-12);
        end

        function test_rotate_90deg_z(testCase)
            q = [cos(pi/4), 0, 0, sin(pi/4)];
            v = [1 0 0];
            v_rot = mr.aux.quat.rotate(q, v);
            testCase.verifyEqual(v_rot, [0 1 0], 'AbsTol', 1e-10);
        end

        function test_rotate_matches_toRotMat(testCase)
            % Verify rotate gives same result as toRotMat * v
            q = mr.aux.quat.normalize([1 2 3 4]);
            v = [5 -3 7];
            v_quat = mr.aux.quat.rotate(q, v);
            R = mr.aux.quat.toRotMat(q);
            v_mat = (R * v.')';
            testCase.verifyEqual(v_quat, v_mat, 'AbsTol', 1e-10);
        end

        function test_rotate_preserves_norm(testCase)
            q = mr.aux.quat.normalize([3 1 4 1]);
            v = [1 2 3];
            v_rot = mr.aux.quat.rotate(q, v);
            testCase.verifyEqual(norm(v_rot), norm(v), 'AbsTol', 1e-10);
        end

    end
end
