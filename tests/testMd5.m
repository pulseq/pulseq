function tests = testMd5
    tests = functiontests(localfunctions);
end

%% Test empty string
function test_empty_string(testCase)
    digest = mr.md5('');
    testCase.verifyEqual(digest, 'd41d8cd98f00b204e9800998ecf8427e');
end

%% Test single character 'a'
function test_single_char_a(testCase)
    digest = mr.md5('a');
    testCase.verifyEqual(digest, '0cc175b9c0f1b6a831c399e269772661');
end

%% Test 'abc'
function test_abc(testCase)
    digest = mr.md5('abc');
    testCase.verifyEqual(digest, '900150983cd24fb0d6963f7d28e17f72');
end

%% Test 'message digest'
function test_message_digest(testCase)
    digest = mr.md5('message digest');
    testCase.verifyEqual(digest, 'f96b697d7cb7938d525a2f31aaf161d0');
end

%% Test 'The quick brown fox jumps over the lazy dog'
function test_quick_brown_fox(testCase)
    digest = mr.md5('The quick brown fox jumps over the lazy dog');
    testCase.verifyEqual(digest, '9e107d9d372bb6826bd81d3542a419d6');
end

%% Test known RFC 1321 vector: alphabet
function test_alphabet(testCase)
    digest = mr.md5('abcdefghijklmnopqrstuvwxyz');
    testCase.verifyEqual(digest, 'c3fcd3d76192e4007dfb496cca67e13b');
end

%% Test determinism (same input → same output)
function test_determinism(testCase)
    d1 = mr.md5('hello world');
    d2 = mr.md5('hello world');
    testCase.verifyEqual(d1, d2);
end
