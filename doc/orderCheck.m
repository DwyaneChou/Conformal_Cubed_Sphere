clc
clear

time_start = 1;
time_end   = 289;

history_path = '..\dynamic_core\run';

res_nc{1} = [history_path,'\','ccsm_output.nc'];

% res_nc{1} = [history_path,'\','ccsm_output_2p0.nc'];
% res_nc{1} = [history_path,'\','ccsm_output_1p0.nc'];
% res_nc{2} = [history_path,'\','ccsm_output_0p5.nc'];
% res_nc{2} = [history_path,'\','ccsm_output.nc'];

res_num = size(res_nc,2);

for ires = 1:res_num
    ics   = ncreadatt(res_nc{ires},'/','ics');
    ice   = ncreadatt(res_nc{ires},'/','ice');
    jcs   = ncreadatt(res_nc{ires},'/','jcs');
    jce   = ncreadatt(res_nc{ires},'/','jce');
    ifs   = ncreadatt(res_nc{ires},'/','ifs');
    ife   = ncreadatt(res_nc{ires},'/','ife');
    xhalo = ncreadatt(res_nc{ires},'/','xhalo');
    yhalo = ncreadatt(res_nc{ires},'/','yhalo');
    
    dx (ires) = ncreadatt(res_nc{ires},'/','dx' );
    dy (ires) = ncreadatt(res_nc{ires},'/','dy' );
    
    its = 1 + xhalo;
    ite = ice;
    jts = 1 + yhalo;
    jte = jce;
    
    areaCell = ncread(res_nc{ires},'areaCell');
    u_end    = ncread(res_nc{ires},'u'      ,[1,1,1,time_end  ],[Inf,Inf,Inf,1]);
    v_end    = ncread(res_nc{ires},'v'      ,[1,1,1,time_end  ],[Inf,Inf,Inf,1]);
    phi_end  = ncread(res_nc{ires},'phi'    ,[1,1,1,time_end  ],[Inf,Inf,Inf,1]);
    u_start  = ncread(res_nc{ires},'u'      ,[1,1,1,time_start],[Inf,Inf,Inf,1]);
    v_start  = ncread(res_nc{ires},'v'      ,[1,1,1,time_start],[Inf,Inf,Inf,1]);
    phi_start= ncread(res_nc{ires},'phi'    ,[1,1,1,time_start],[Inf,Inf,Inf,1]);
    
    L1_phi  (ires) = L1  (phi_end(its:ite,jts:jte,:),phi_start(its:ite,jts:jte,:),areaCell(its:ite,jts:jte,:));
    L2_phi  (ires) = L2  (phi_end(its:ite,jts:jte,:),phi_start(its:ite,jts:jte,:),areaCell(its:ite,jts:jte,:));
    LInf_phi(ires) = LInf(phi_end(its:ite,jts:jte,:),phi_start(its:ite,jts:jte,:),areaCell(its:ite,jts:jte,:));
end

for ires = 2:res_num
    L1order_phi  (ires) = log(L1_phi  (ires)/L1_phi  (ires-1))/log(dx(ires)/dx(ires-1));
    L2order_phi  (ires) = log(L2_phi  (ires)/L2_phi  (ires-1))/log(dx(ires)/dx(ires-1));
    LInforder_phi(ires) = log(LInf_phi(ires)/LInf_phi(ires-1))/log(dx(ires)/dx(ires-1));
end

function reslut = L1(field_model,field_ref,areaCell)

reslut = sum(sum(sum(abs(field_model - field_ref) .* areaCell)))...
       / sum(sum(sum(abs(field_ref) .* areaCell)));

end

function reslut = L2(field_model,field_ref,areaCell)

reslut = sqrt(sum(sum(sum((field_model - field_ref).^2 .* areaCell)))...
             /sum(sum(sum(field_ref.^2 .* areaCell))));

end

function reslut = LInf(field_model,field_ref,areaCell)

reslut = max(max(max(abs(field_model - field_ref) .* areaCell)))...
       / max(max(max(abs(field_ref) .* areaCell)));

end