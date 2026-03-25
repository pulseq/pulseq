classdef InputParserCompat < inputParser
  % this is a customized parser class to replicate Matlab's undocumented
  % behaviour on GNU Octave. In Matlab the optional positional arguments can
  % also be called as param-value-pair with an arbitrary position
  properties
    needCopyOptional = mr.aux.isOctave();
  end

  methods
    function parse (this, varargin)
      if this.needCopyOptional && length(this.Optional)>0
        for i=1:length(this.Optional)
          this.Parameter(end+1)=this.Optional{i};
        end
        this.needCopyOptional = false;
      end
      parse@inputParser(this, varargin{:});
    end
  end
end
