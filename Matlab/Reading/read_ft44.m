function [neut,wld] = read_ft44(file,quiet)
% [neut,wld] = read_ft44(file)
%
% Read fort.44 file
% 
% J.D. Lore, modified from W. Dekeyser routine

if nargin < 2
    quiet = 0;
end

nlwrmsh = 1;  

[fid,msg] = fopen(file);
if (fid == -1)
    error(msg);
end

% neut = struct([]);
wld  = struct;

%% Read dimensions

% nx, ny, version
dims = fscanf(fid,'%d',3);
nx   = dims(1);
ny   = dims(2);
ver  = dims(3);

if ver ~= 20081111 && ver ~= 20160829 && ver ~= 20170328 && ~quiet
    %     fprintf('Warning: read_ft44: unknown format of fort.44 file :%d\n',ver);
end

% go to new line (skip reading a possible git-hash)
fgetl(fid);

% natm, nmol, nion
dims = fscanf(fid,'%d',3);
natm = dims(1);
nmol = dims(2);
nion = dims(3);

neut.natm = natm;
neut.nmol = nmol;
neut.nion = nion;

line = fgetl(fid);
for i = 1:natm
    neut.atm{i} = strtrim(fgetl(fid));
end
for i = 1:nmol
    neut.mol{i} = strtrim(fgetl(fid));
end
for i = 1:nion
    neut.ion{i} = strtrim(fgetl(fid));
end
neut.speciesArray = [neut.atm,neut.mol,neut.ion];

%% Read basic data

neut.dab2     = read_ft44_rfield(fid,ver,'dab2',[nx,ny,natm]);
neut.tab2     = read_ft44_rfield(fid,ver,'tab2',[nx,ny,natm]);
neut.dmb2     = read_ft44_rfield(fid,ver,'dmb2',[nx,ny,nmol]);
neut.tmb2     = read_ft44_rfield(fid,ver,'tmb2',[nx,ny,nmol]);
neut.dib2     = read_ft44_rfield(fid,ver,'dib2',[nx,ny,nion]);
neut.tib2     = read_ft44_rfield(fid,ver,'tib2',[nx,ny,nion]);
neut.rfluxa   = read_ft44_rfield(fid,ver,'rfluxa',[nx,ny,natm]);
neut.rfluxm   = read_ft44_rfield(fid,ver,'rfluxm',[nx,ny,nmol]);
neut.pfluxa   = read_ft44_rfield(fid,ver,'pfluxa',[nx,ny,natm]);
neut.pfluxm   = read_ft44_rfield(fid,ver,'pfluxm',[nx,ny,nmol]);
neut.refluxa  = read_ft44_rfield(fid,ver,'refluxa',[nx,ny,natm]);
neut.refluxm  = read_ft44_rfield(fid,ver,'refluxm',[nx,ny,nmol]);
neut.pefluxa  = read_ft44_rfield(fid,ver,'pefluxa',[nx,ny,natm]);
neut.pefluxm  = read_ft44_rfield(fid,ver,'pefluxm',[nx,ny,nmol]);
neut.emiss    = read_ft44_rfield(fid,ver,'emiss',[nx,ny,1]);
neut.emissmol = read_ft44_rfield(fid,ver,'emissmol',[nx,ny,1]);
neut.srcml    = read_ft44_rfield(fid,ver,'srcml',[nx,ny,nmol]);
neut.edissml  = read_ft44_rfield(fid,ver,'edissml',[nx,ny,nmol]);
% % % %% Data on wall loading
% % %
% nlim, nsts, nstra
% warning('need a check if no wlld data')
dims  = fscanf(fid,'%d',3);
nlim  = dims(1);
nsts  = dims(2);
nstra = dims(3);
wld.nlim = nlim;
wld.nsts = nsts;
wld.nstra = nstra;
wld.wldnek = zeros(nlim+nsts,nstra+1);
wld.wldnep = zeros(nlim+nsts,nstra+1);
wld.wldna  = zeros(nlim+nsts,natm,nstra+1);
wld.ewlda  = zeros(nlim+nsts,natm,nstra+1);
wld.wldnm  = zeros(nlim+nsts,nmol,nstra+1);
wld.ewldm  = zeros(nlim+nsts,nmol,nstra+1);
% % %
% try
wld.wldnek(:,1)  = read_ft44_rfield(fid,ver,'wldnek',[nlim+nsts]);
wld.wldnep(:,1)  = read_ft44_rfield(fid,ver,'wldnep',nlim+nsts);
wld.wldna(:,:,1) = read_ft44_rfield(fid,ver,'wldna',[nlim+nsts,natm]);
wld.ewlda(:,:,1) = read_ft44_rfield(fid,ver,'ewlda',[nlim+nsts,natm]);
wld.wldnm(:,:,1) = read_ft44_rfield(fid,ver,'wldnm',[nlim+nsts,nmol]);
wld.ewldm(:,:,1) = read_ft44_rfield(fid,ver,'ewldm',[nlim+nsts,nmol]);
% end
% wall_Geometry
wld.poly = read_ft44_rfield(fid,ver,'wall_geometry',[4,nlim]);

wld.wldra  = zeros(nlim+nsts,natm,nstra+1);
wld.wldrm  = zeros(nlim+nsts,nmol,nstra+1);

wld.wldra(:,:,1)  = read_ft44_rfield(fid,ver,'wldra',[nlim+nsts,natm]);
wld.wldrm(:,:,1)  = read_ft44_rfield(fid,ver,'wldrm',[nlim+nsts,nmol]);

if (nstra > 1)
    for i = 2:nstra+1
        wld.wldnek(:,i)  = read_ft44_rfield(fid,ver,'wldnek',nlim+nsts);
        wld.wldnep(:,i)  = read_ft44_rfield(fid,ver,'wldnep',nlim+nsts);
        wld.wldna(:,:,i) = read_ft44_rfield(fid,ver,'wldna',[nlim+nsts,natm]);
        wld.ewlda(:,:,i) = read_ft44_rfield(fid,ver,'ewlda',[nlim+nsts,natm]);
        wld.wldnm(:,:,i) = read_ft44_rfield(fid,ver,'wldnm',[nlim+nsts,nmol]);
        wld.ewldm(:,:,i) = read_ft44_rfield(fid,ver,'ewldm',[nlim+nsts,nmol]);
        wld.wldra(:,:,i)  = read_ft44_rfield(fid,ver,'wldra',[nlim+nsts,natm]);
        wld.wldrm(:,:,i)  = read_ft44_rfield(fid,ver,'wldrm',[nlim+nsts,nmol]);
    end
end

% Need to figure out NLPS
[wld.wldpp(:,:,1),dims_out] = read_ft44_rfield(fid,ver,'wldpp',[nlim+nsts,1],[0,1]);

npls = dims_out(2);
wld.npls = npls;

wld.wldpp   = zeros(nlim+nsts,npls,nstra+1);
wld.wldpa   = zeros(nlim+nsts,natm,nstra+1);
wld.wldpm   = zeros(nlim+nsts,nmol,nstra+1);
wld.wldpeb  = zeros(nlim+nsts,nstra+1);
wld.wldspt  = zeros(nlim+nsts,nstra+1);
wld.wldspta  = zeros(nlim+nsts,natm,nstra+1);
wld.wldsptm  = zeros(nlim+nsts,nmol,nstra+1);


wld.wldpa(:,:,1) = read_ft44_rfield(fid,ver,'wldpa',[nlim+nsts,natm]);
wld.wldpm(:,:,1) = read_ft44_rfield(fid,ver,'wldpm',[nlim+nsts,nmol]);
wld.wldpeb(:,1)  = read_ft44_rfield(fid,ver,'wldpeb',nlim+nsts);
wld.wldspt(:,1)  = read_ft44_rfield(fid,ver,'wldspt',nlim+nsts);
if ver > 20080706
    wld.wldspta(:,:,1) = read_ft44_rfield(fid,ver,'wldspta',[nlim+nsts,natm]);
    wld.wldsptm(:,:,1)  = read_ft44_rfield(fid,ver,'wldsptm',[nlim+nsts,nmol]);
end
if (nstra > 1)
    for i = 2:nstra+1
        wld.wldpp(:,:,i) = read_ft44_rfield(fid,ver,'wldpp',[nlim+nsts,npls]);
        wld.wldpa(:,:,i) = read_ft44_rfield(fid,ver,'wldpa',[nlim+nsts,natm]);
        wld.wldpm(:,:,i) = read_ft44_rfield(fid,ver,'wldpm',[nlim+nsts,nmol]);
        wld.wldpeb(:,i)  = read_ft44_rfield(fid,ver,'wldpeb',nlim+nsts);
        wld.wldspt(:,i)  = read_ft44_rfield(fid,ver,'wldspt',nlim+nsts);
        if ver > 20080706
            wld.wldspta(:,:,i) = read_ft44_rfield(fid,ver,'wldspta',[nlim+nsts,natm]);
            wld.wldsptm(:,:,i)  = read_ft44_rfield(fid,ver,'wldsptm',[nlim+nsts,nmol]);
        end
    end
end

wld.isrftype = read_ft44_rfield(fid,ver,'isrftype',nlim+nsts);


if nlwrmsh
    
    wld.wlarea  = read_ft44_rfield(fid,ver,'wlarea',nlim+nsts,[],1);
    
    wld.wlabsrp = zeros(natm+nmol+nion+npls,nlim+nsts);
    wld.wlabsrp(1:natm,:) = read_ft44_rfield(fid,ver,'wlabsrp',[natm,nlim+nsts],[],ceil(natm/6));
    ns = natm;
    wld.wlabsrp(ns+1:ns+nmol,:) = read_ft44_rfield(fid,ver,'wlabsrp',[nmol,nlim+nsts],[],ceil(nmol/6));
    ns = natm+nmol;
    wld.wlabsrp(ns+1:ns+nion,:) = read_ft44_rfield(fid,ver,'wlabsrp',[nion,nlim+nsts],[],ceil(nion/6));
    ns = natm+nmol+nion;
    wld.wlabsrp(ns+1:ns+npls,:) = read_ft44_rfield(fid,ver,'wlabsrp',[npls,nlim+nsts],[],ceil(npls/6));
    
    % take this logic to wlabsrp too!
    
    wld.wlpump = zeros(natm+nmol+nion+npls,nlim+nsts);
    if ver > 20080706
        istart = 1;
        ns = natm;
        if ns > 0
            wld.wlpump(istart:istart+ns-1,:) = read_ft44_rfield(fid,ver,'wlpump',[ns,nlim+nsts],[],ceil(ns/6));
            istart = istart + ns;
        end
        ns = nmol;
        if ns > 0
            wld.wlpump(istart:istart+ns-1,:) = read_ft44_rfield(fid,ver,'wlpump',[ns,nlim+nsts],[],ceil(ns/6));
            istart = istart + ns;
        end
        ns = nion;
        if ns > 0
            wld.wlpump(istart:istart+ns-1,:) = read_ft44_rfield(fid,ver,'wlpump',[ns,nlim+nsts],[],ceil(ns/6));
            istart = istart + ns;
        end
        ns = npls;
        if ns > 0
            wld.wlpump(istart:istart+ns-1,:) = read_ft44_rfield(fid,ver,'wlpump',[ns,nlim+nsts],[],ceil(ns/6));
        end
    end
    if natm > 0
        neut.eneutrad = read_ft44_rfield(fid,ver,'eneutrad',[nx,ny,natm]);
    end
    if ver > 20080706
        if nmol > 0
            neut.emolrad = read_ft44_rfield(fid,ver,'emolrad',[nx,ny,nmol]);
        end
        if nion > 0
            neut.eionrad = read_ft44_rfield(fid,ver,'eionrad',[nx,ny,nion]);
        end
    end
    
    % C Meaning of the arrays which are written to the fort.44 file
    % C
    % C  eirdiag_nds_typ    : type of the surface
    % C                       1: poloidal surface (Y=CONST)
    % C                       2: radial surface (X=CONST)
    % C  eirdiag_nds_srf    : index of the corresponding radial or poloidal SURFACE
    % C                       on which the distribution is given
    % C  eirdiag_nds_start  : index of the first CELL on the surface '..._srf'
    % C  eirdiag_nds_end    : ... of the last CELL ...
    % C
    % C  _sfr, _start and _end are given in B2 indexing.
    % C  '_srf' corresponds to the indexing of B2 arrays with fluxes (CELL-1 for outer wall and target).
    % C
    % C  eirdiag_nds_ind(IS) : is the index AFTER WHICH the data for surface IS is started in arrays '*_res'.
    % C                        IF _ind(IS)<0 then this surface is skipped
    % C
    
    wld.eirdiag = read_ft44_rfield(fid,ver,'eirdiag',5*nsts+1);
    istart = 1; iend = nsts + 1;
    wld.eirdiag_nds_ind = wld.eirdiag(istart:iend);
    istart = iend + 1; iend = istart + nsts-1;
    wld.eirdiag_nds_typ  = wld.eirdiag(istart:iend);
    istart = iend + 1; iend = istart + nsts-1;
    wld.eirdiag_nds_srf  = wld.eirdiag(istart:iend);
    istart = iend + 1; iend = istart + nsts-1;
    wld.eirdiag_nds_start  = wld.eirdiag(istart:iend);
    istart = iend + 1; iend = istart + nsts-1;
    wld.eirdiag_nds_end    = wld.eirdiag(istart:iend);
    
    ncl = wld.eirdiag_nds_ind(nsts+1);
    wld.ncl = ncl;
    % % I think these go [total,section,total,section,...]    
    wld.sarea_res  = read_ft44_rfield(fid,ver,'sarea_res',ncl);
    wld.wldna_res  = read_ft44_rfield(fid,ver,'wldna_res',[natm,ncl]);
    wld.wldnm_res  = read_ft44_rfield(fid,ver,'wldnm_res',[nmol,ncl]);
    wld.ewlda_res = read_ft44_rfield(fid,ver,'ewlda_res',[natm,ncl]);
    wld.ewldm_res = read_ft44_rfield(fid,ver,'ewldm_res',[nmol,ncl]);
    wld.ewldea_res = read_ft44_rfield(fid,ver,'ewldea_res',[natm,ncl]);
    wld.ewldem_res = read_ft44_rfield(fid,ver,'ewldem_res',[nmol,ncl]);
    wld.ewldrp_res = read_ft44_rfield(fid,ver,'ewldrp_res',ncl);
    wld.ewldmr_res = read_ft44_rfield(fid,ver,'ewldmr_res',[nmol,ncl]);
    
    if ver > 20080706
        wld.wldspt_res  = read_ft44_rfield(fid,ver,'wldspt_res',ncl);
        wld.wldspta_res  = read_ft44_rfield(fid,ver,'wldspta_res',[ncl,natm]);
        wld.wldsptm_res  = read_ft44_rfield(fid,ver,'wldsptm_res',[ncl,nmol]);
        
        istart = 1;
        wld.wlpump_res = zeros(ncl,natm+nmol+nion+npls);
        ns = natm;
        if ns > 0
            wld.wlpump_res(:,istart:istart+ns-1) = read_ft44_rfield(fid,ver,'wlpump_res',[ncl,ns],[],ceil(ns/6));
            istart = istart + ns;
        end
        ns = nmol;
        if ns > 0
            wld.wlpump_res(:,istart:istart+ns-1) = read_ft44_rfield(fid,ver,'wlpump_res',[ncl,ns],[],ceil(ns/6));
            istart = istart + ns;
        end
        ns = nion;
        if ns > 0
            wld.wlpump_res(:,istart:istart+ns-1) = read_ft44_rfield(fid,ver,'wlpump_res',[ncl,ns],[],ceil(ns/6));
            istart = istart + ns;
        end
        ns = npls;
        if ns > 0
            wld.wlpump_res(:,istart:istart+ns-1) = read_ft44_rfield(fid,ver,'wlpump_res',[ncl,ns],[],ceil(ns/6));
        end
        
        
        % %
        neut.pdena_int    = read_ft44_rfield(fid,ver,'pdena_int',[natm,nstra+1]);
        neut.pdenm_int    = read_ft44_rfield(fid,ver,'pdenm_int',[nmol,nstra+1]);
        neut.pdeni_int    = read_ft44_rfield(fid,ver,'pdeni_int',[nion,nstra+1]);
        neut.pdena_int_b2 = read_ft44_rfield(fid,ver,'pdena_int_b2',[natm,nstra+1]);
        neut.pdenm_int_b2 = read_ft44_rfield(fid,ver,'pdenm_int_b2',[nmol,nstra+1]);
        neut.pdeni_int_b2 = read_ft44_rfield(fid,ver,'pdeni_int_b2',[nion,nstra+1]);
        neut.edena_int    = read_ft44_rfield(fid,ver,'edena_int',[natm,nstra+1]);
        neut.edenm_int    = read_ft44_rfield(fid,ver,'edenm_int',[nmol,nstra+1]);
        neut.edeni_int    = read_ft44_rfield(fid,ver,'edeni_int',[nion,nstra+1]);
        neut.edena_int_b2 = read_ft44_rfield(fid,ver,'edena_int_b2',[natm,nstra+1]);
        neut.edenm_int_b2 = read_ft44_rfield(fid,ver,'edenm_int_b2',[nmol,nstra+1]);
        neut.edeni_int_b2 = read_ft44_rfield(fid,ver,'edeni_int_b2',[nion,nstra+1]);
        % % else
        % %     neut.eneutrad = read_ft44_rfield(fid,[nx,ny,natm]);
    end
end



% figure; hold on;
% for i = 1:nlim
%     istart = (i-1)*4 + 1;
%     plot([wld.poly(istart),wld.poly(istart+2)],[wld.poly(istart+1),wld.poly(istart+3)])
% end




fclose(fid);
