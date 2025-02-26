function OUT = isOctave ()
% function that checks if we are in Octave
  persistent IO;
  if (isempty (IO))
    IO = exist ('OCTAVE_VERSION', 'builtin');
  end
  OUT = IO;
end
