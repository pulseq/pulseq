function [melody, timeSign] = melodyToScale(melody,speed)
% mrMusic.melodyToPitchesAdnDurations : convert *melody* to another
% *melody* object that consists of all notes used in the original melody
% played one after one another in the ascending order to form a scale-like
% pattern. It is mainly useful for checking the resonances of the
% particular scanner and making sure that all notes sound equaly nice on
% the real hardware. As a duration of the note the typical note duration
% from the melody is used. The optional paramater *speed* allows to slow 
% down (speed < 1) or accelerate (speed>1) the tempo.

allnotes=[melody{:}];
alldurs=imag(allnotes);
td=median(alldurs); % typical duration

if nargin>1
    td=td/speed;
end

allfreqs=unique(real(allnotes./alldurs));
allfreqs=allfreqs(allfreqs>0);
melody={...
    td*[allfreqs+1i],td*[zeros(size(allfreqs))+1i],td*[zeros(size(allfreqs))+1i];...
    td*[zeros(size(allfreqs))+1i],td*[allfreqs+1i],td*[zeros(size(allfreqs))+1i];...
    td*[zeros(size(allfreqs))+1i],td*[zeros(size(allfreqs))+1i],td*[allfreqs+1i];...
    };
timeSign=td*length(allfreqs);
