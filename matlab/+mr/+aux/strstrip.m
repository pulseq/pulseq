function out=strstrip(in)
%STRSTRIP: clean up the beginning and the end of a character string
  if ~mr.aux.isOctave
    if ~isempty(in) 
        out=strip(in);
    else
        out=in;
    end
  else
    % octave doen't have strip
    while ~isempty(in) && (in(1)==' ' || in(1)=="\t" || in(1)=="\r" || in(1)=="\n")
      in(1)=[];
    end
    while ~isempty(in) && (in(end)==' ' || in(end)=="\t" || in(end)=="\r" || in(end)=="\n")
      in(end)=[];
    end
    out=in;
  end
end
