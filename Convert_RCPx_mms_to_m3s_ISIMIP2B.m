%% Convert ISIMIP2B runoff from kg m-2 s (mm/day) to m3s
clear all

root = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
root_data = sprintf('%s\\data\\data_prep\\ISIMIP_Runoff_ISIMIP2B\\mms',root);

%% Data list
fnames = dir(root_data);

for j=1:numel(fnames)
    
    R_names{j} = sprintf('%s\\%s',root_data,(fnames(j).name));
end

%% lat-map
fname = sprintf('%s\\data\\data_prep\\Latmap\\latmap.mat',root);
load(fname)

%%
for k=1:numel(R_names)
    clear R Rm3s R_data
    if k==1 || k==2; continue; end;
    
    fprintf('Reading file %d of #%d\n',k,numel(R_names));
    
    R_data = ncread(R_names{k},'qtot');
    
    %% Data prep
    disp('Prep')
    mv = linspace(11,0,12);
    
    for m=1:12
        R{m} = R_data(:,:,m)';
    end
    
    %% Yearly average
    disp('Yearly average')
    [nr,nc]=size(R{1});
    for r=1:nr
        for c=1:nc
            R{13}(r,c) = (R{1}(r,c)+R{2}(r,c)+R{3}(r,c)+R{4}(r,c)+R{5}(r,c)+R{6}(r,c)+R{7}(r,c)+R{8}(r,c)+R{9}(r,c)+R{10}(r,c)+R{11}(r,c)+R{12}(r,c))/12;
        end
    end
    
    %% Convert kg m-2 s-1 (mm/s) to m3/s
    disp('Converting to m3s')
    for m=1:13
        R{m}(find(isnan(R{m}(:))))=0;
        
        for r=1:nr
            for c=1:nc
                Rm3s{m}(r,c) = R{m}(r,c) * 1e-3 * (111e3*0.5)^2 * cos(degtorad(latmap(r,c)));
            end
        end
        
    end
    
    %% Save 30yr average m3s map
    disp('Saving')
    matfile = fullfile(root, sprintf('data\\data_prep\\ISIMIP_Runoff_ISIMIP2B\\m3s\\%s_m3s.mat',fnames(k).name(1:(end-4))));
    save(matfile,'Rm3s','-v7.3');
end