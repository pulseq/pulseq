classdef EventLibrary < handle
    % EventLibrary   Maintain a list of events.
    %
    % The class is used by the Sequence class to store events of an MRI 
    % sequence defined using the Pulseq file format.
    %   See http://pulseq.github.io/
    %
    % This class is designed for *performance* so some functions can be
    % compiled to mex files using MATLAB's code generation. A parallel list
    % of IDs, data, and meta information is maintained to maximise search
    % speed.
    %
    % Sequence Properties:
    %    keys - A list of event IDs
    %    data - A struct array with field 'array' to store data of varying
    %           lengths, remaining compatible with codegen.
    %    lengths - Corresponding lengths of the data arrays
    %    type - Type to distinguish events in the same class (e.g.
    %           trapezoids and arbitrary gradients)
    %
    % Sequence Methods:
    %    find - Find an event in the library
    %    insert - Add a new event to the library
    %
    % See also   mr.Sequence
    %
    % Kelvin Layton <kelvin.layton@uniklinik-freiburg.de>
    
    properties
        keys;
        data;
        lengths;
        type;
        keymap;
    end
    
    methods
        function obj = EventLibrary()
            obj.keys=zeros(1,0);
            obj.data=struct('array',{});
            obj.lengths=zeros(1,0);
            obj.type=char(zeros(1,0));
            obj.keymap=containers.Map('KeyType', 'char', 'ValueType', 'double'); % unsure if int64 is a good choice
        end
        
        function [id, found] = find(obj,data)
            %find Lookup a data structure in the given library.
            %   idx=find(lib,data) Return the index of the data in
            %   the library. If the data does not exist in the library then
            %   the index for the next new entry is returned.
            %
            %   The data is a 1xN array with event-specific data.
            %
            %   See also  insert mr.Sequence.addBlock
            
            %try
            %    [id, found] = mr.EventLibrary.find_mex(obj.keys,obj.data,obj.lengths,data);
            %catch 
            %    [id, found] = mr.EventLibrary.find_mat(obj.keys,obj.data,obj.lengths,obj.keymap,data);
            %end
            
            % use map index for faster searches
            % matlab is extremely limited with regard to advanced contasiners
            % we therefore are forced to use hashed map and convert data to a
            % string
            data_string=sprintf('%.6g ', data); % precision can be discussed
            % containers.Map does not have a proper find function so we use direct
            % access and catch the possible error
            try
                id=obj.keymap(data_string);
                found=true;
            catch 
                if isempty(obj.keys)
                    id=1;
                else
                    id=max(obj.keys)+1;
                end
                found=false;
            end
        end
               
        
        function insert(obj,id,data,type)
            %insert Add event to library
            % 
            % See also find
            
            obj.keys(id)=id;
            obj.data(id).array=data;
            obj.lengths(id)=length(data);
            
            % use map index for faster searches
            % matlab is extremely limited with regard to advanced contasiners
            % we therefore are forced to use hashed map and convert data to a
            % string
            %data_string=num2hex(data)';
            %data_string=data_string(:)';
            data_string=sprintf('%.6g ', data);
            obj.keymap(data_string)=id;
            
            if nargin>3
                obj.type(id)=type;
            end
        end
        
    end
    
    methods(Static)
        % Helper functions for fast searching
        % See find_mat.m
        [id,found] = find_mat(keys,data,lengths,keymap,newData);
        %[id,found] = find_mex(keys,data,lengths,newData);
    end
    
end

