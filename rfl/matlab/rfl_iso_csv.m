function [R,E_GRID,R_UNIT,E_UNIT,inst_info] = rfl_iso_csv(inst_info,species,outfile,channels,varargin)
% [R,E_GRID,R_UNIT,E_UNIT,inst_info] = rfl_iso_csv(inst_info,species,outfile,channels,...)
% write isotropic response to csv output file
% inst_info  - string giving file name or structure containing inst_info
% species - string giving species to study (e.g., 'PROT')
% outfile - output file name (absent or empty to not write file)
% channels - cell array of channel names (absent or empty to do all channels)
% options: (provided as name,value pairs at end of function call)
%  'comments',true - include extrea headers like species and inst_info REFERENCES
%  'comment_char','#' - preceed extra headers with a special character, like '#'
% OUTPUTS:
% R - NE x Nchannels isotropic response
% E_GRID - energy
% R_UNIT - units for R (usually cm^2 sr)
% E_UNIT - units of E_GRID (usually MeV)
% inst_info - loaded inst_info structure

inst_info = rfl_load_inst_info(inst_info);

if ~ismember(species,inst_info.SPECIES)
    error('Species %s not found in inst_info',species);
end

if nargin < 3
    outfile = '';
end

if nargin < 4
    channels = {};
end

comments = false;
comment_char = '';

for i = 1:2:length(varargin)
    switch(lower(varargin{i}))
        case 'comments'
            comments = varargin{i+1};
        case 'comment_char'
            comment_char = varargin{i+1};
        otherwise
            error('Unknown argument %s',varargin{i});
    end
end

if isempty(channels)
    channels = inst_info.CHANNEL_NAMES;
end

R = [];
E_GRID = [];
E_UNIT = inst_info.E_UNIT;
L_UNIT = inst_info.L_UNIT;
kept_channels = {};
for ichan = 1:length(channels)
    chan = channels{ichan};
    if ~ismember(species,inst_info.(chan).SPECIES)
        warning('Skipping %s, has no entry for %s',chan,species);
        continue; % skip channel, doesn't have response for requested species
    end
    resp = inst_info.(chan).(species);
    if ~isequal(resp.E_UNIT,E_UNIT)
        error('Energy unit %s for %s/%s inconsistent with main file, %s',resp.L_UNIT,chan,species,L_UNIT);
    end
    if ~isequal(resp.L_UNIT,L_UNIT)
        error('Length unit %s for %s/%s inconsistent with main file, %s',resp.L_UNIT,chan,species,L_UNIT);
    end
    if isempty(E_GRID)
        E_GRID = resp.E_GRID;
    else
        if (~isequal(size(E_GRID),size(resp.E_GRID))) || (max(abs(resp.E_GRID./E_GRID-1))>0.0001)
            error('Energy grids are not the same across all requested channels (devaition found in %s/%s)\n',chan,species);
        end
    end
    
    hE = resp.make_hE(resp,resp.E_GRID,[]); % hE integrates over angle but not energy
    Ri = hE./rfl_make_deltas(resp.E_GRID); % remove delta-E to just get the isotropic energy response vector
    kept_channels{end+1} = chan;
    R = [R,Ri];
end

R_UNIT = sprintf('%s^2 sr',L_UNIT);

if ~isempty(outfile) % write to CSV
    if isempty(kept_channels)
        warning('No requested channels had response for %s. Not writing to %s',species,outfile);
    else
        fid = -1;
        try
            fid = fopen(outfile,'w');
            if comments
                fprintf(fid,'%sIsotropic response file produced by Response Function Library\n',comment_char);
                fprintf(fid,'%sSPECIES: %s\n',comment_char,species);
                fprintf(fid,'%sR unit: %s\n',comment_char,R_UNIT);
                fprintf(fid,'%sBased on:\n',comment_char);
                for i = 1:length(inst_info.REFERENCES)
                    fprintf(fid,'%s%s\n',comment_char,inst_info.REFERENCES{i});
                end
            end
            fprintf(fid,'Energy %s',E_UNIT);
            fprintf(fid,',%s',kept_channels{:});
            if ~comments
                fprintf(fid,' %s',R_UNIT); % append R_UNIT to last channel header
            end
            fprintf(fid,'\n');
            out = [E_GRID(:),R];
            fmt = repmat('%g,',1,size(out,2));
            fmt(end+[0:1]) = '\n';
            fprintf(fid,fmt,out');
            fclose(fid);
        catch ME
            if ~isequal(fid,-1)
                fclose(fid);
            end
            rethrow(ME)
        end
    end
end
