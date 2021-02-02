clearvars;

fullPathToRun = '/Example';

Geo   = read_b2fgmtry_simple(fullfile(fullPathToRun,'b2fgmtry'));
State = read_b2fstate_simple(fullfile(fullPathToRun,'b2fstate'));

cellR = reshape(Geo.crx(:,:,[1,2,4,3]),(Geo.nx+2)*(Geo.ny+2),4);
cellZ = reshape(Geo.cry(:,:,[1,2,4,3]),(Geo.nx+2)*(Geo.ny+2),4);

figure;
patch(cellR.',cellZ.',State.te(:),'edgecolor','none')
