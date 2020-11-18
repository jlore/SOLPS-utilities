function Geo = read_b2fgmtry_simple(file,quiet)
% Geo = read_b2fgmtry_simple(file,quiet)
%
% Reads b2fgmtry, and infers dimensions
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

%% Get version of the b2fgmtry file
line    = fgetl(fid);
version = line(8:17);
if ~quiet
    disp(['read_b2fgmtry_simple -- file version ',version]);
end
Geo.version = version;

%% Get dimensions and sanity check
line    = fgetl(fid);
data = textscan(line,'%s %s %d %s');
assert(strcmp(data{2},'int') && data{3}==2 && strcmp(data{4},'nx,ny'),'Do not understand format of b2fgmtry')
data = fscanf(fid,'%d',2);
Geo.nx = data(1);
Geo.ny = data(2);

numCells = (Geo.nx+2)*(Geo.ny+2);
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
        
        % Check if this is a (nx+2,ny+2,dim3) sized array
        if mod(varSize,numCells) ==0
            this = reshape(this,[Geo.nx+2,Geo.ny+2,varSize/numCells]);
        end
        Geo.(varName) = this;
    end
end

fclose(fid);



