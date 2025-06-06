function [version_major, version_minor, version_revision, version_combined]=version(type)
%version: return the relevant version information. The type of the version
%         information is defined by the (optional) parameter.
%         Possible version types: 
%    'pulseq' : (default) verson of the current Matlab package
%    'output' : version of the file written by the seq.write() function

if nargin==0 || strcmp(type,'pulseq')
    version_major=1;
    version_minor=5;
    version_revision=1;
elseif strcmp(type,'output')
    version_major=1;
    version_minor=5;
    version_revision=1;
else
    error('Unsupported version request, type=%s',type);
end
version_combined=1000000*version_major+1000*version_minor+version_revision;
