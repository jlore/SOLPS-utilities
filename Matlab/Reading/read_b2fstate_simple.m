function State = read_b2fstate_simple(file,quiet)
% State = read_b2fstate_simple(file,quiet)
%
% Reads b2fstat[e,i], and infers dimensions
%
% J.D. Lore, 2020

if nargin < 2
    quiet = 0;
end

%% Open file
[fid,msg] = fopen(file);
if (fid == -1)
    error(msg);
end

%% Get version of the state file
line    = fgetl(fid);
version = line(8:17);
if ~quiet
    disp(['read_b2fstate_simple -- file version ',version]);
end
State.version = version;

%% Get dimensions and sanity check
line    = fgetl(fid);
data = textscan(line,'%s %s %d %s');
assert(strcmp(data{2},'int') && data{3}==3 && strcmp(data{4},'nx,ny,ns'),'Do not understand format of b2fstate')
data = fscanf(fid,'%d',3);
State.nx = data(1);
State.ny = data(2);
State.ns = data(3);

if State.ns == 1
    warning("Not sure this will work for a single species!")
end

numCells = (State.nx+2)*(State.ny+2);
while ~feof(fid)
    
    line=fgets(fid);
    
    if strcmp(line(1),'*')
        data = textscan(line,'%s %s %d %s');
        varType = char(data{2});
        varSize = data{3};
        varName = char(data{4});
        switch varType
            case 'char'
                this = fgetl(fid);
            case 'int'
                this = fscanf(fid,'%d',varSize);
            case 'real'
                this = fscanf(fid,'%e',varSize);
            otherwise
                error('Do not understand data format: %s',varType)
        end

        % Test for size of array and try to deduce dimensions
        switch varSize
            case numCells
                % This is a scalar quantity
                this = reshape(this,[State.nx+2,State.ny+2]);
            case numCells*2
                % This is a flux quantity
                this = reshape(this,[State.nx+2,State.ny+2,2]);
            case numCells*State.ns
                % This is a scalar quantity by species
                this = reshape(this,[State.nx+2,State.ny+2,State.ns]);
            case numCells*State.ns*2
                % This is a flux quantity by species
                this = reshape(this,[State.nx+2,State.ny+2,2,State.ns]);
            case numCells*4
                % This is a flux quantity by species in 3.1 format
                this = reshape(this,[State.nx+2,State.ny+2,2,2]);
            case numCells*State.ns*4                
                % This is a flux quantity by species in 3.1 format
                this = reshape(this,[State.nx+2,State.ny+2,2,2,State.ns]);                
            otherwise
                if mod(varSize,numCells) == 0
                    warning('Must have missed some dimension checks above for line %s',line)
                end
        end
        
        State.(varName) = this;
    end
end

fclose(fid);



