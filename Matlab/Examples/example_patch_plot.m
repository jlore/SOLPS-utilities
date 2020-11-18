clearvars;

Geo   = read_b2fgmtry_simple('C:/Users/jjl/Dropbox (ORNL)/SPARC_SOLPS/V1E_LSN2_D+C/baserun/b2fgmtry');
State = read_b2fstate_simple('C:/Users/jjl/Dropbox (ORNL)/SPARC_SOLPS/V1E_LSN2_D+C/P29MW_n1.65e20_Rout0.9_Y2pc/b2fstate');

cellR = reshape(Geo.crx(:,:,[1,2,4,3]),(Geo.nx+2)*(Geo.ny+2),4);
cellZ = reshape(Geo.cry(:,:,[1,2,4,3]),(Geo.nx+2)*(Geo.ny+2),4);

figure;
patch(cellR.',cellZ.',State.te(:),'edgecolor','none')