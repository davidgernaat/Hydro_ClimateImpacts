% downscale LPJ runoff R to 15s and 3s resolutions
clear all;
res=15;
root = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
CONT={'AFR','ASIA','AUS','CAM','EUR','NAM','SAM'};
nc=[1 2 3 4 5 6 7];

%% File names
% for ir=1
ir=4;
clearvars -except root RCP_name ir res CONT nc

root_data = sprintf('%s\\data\\data_prep\\ISIMIP_Runoff_ISIMIP2B\\m3s',root);

%% Data list
fnames = dir(root_data);

for j=1:numel(fnames)
    R_names{j} = sprintf('%s\\%s',root_data,(fnames(j).name));
end

%%
for y=1:numel(R_names)
    if y==1 || y==2; continue; end;
    fprintf('\nReading %d of #%d \n',y, numel(R_names))
    
    load(R_names{y});
    R_data = Rm3s;
    
    %%
    for nt=nc
        clearvars -except res root CONT nc acc y R_data nt R_names R_name ir fnames RCP_name
        
        fprintf('Continent: %s\n',CONT{nt})
        
        %disp('Loading ACC');
        fname = sprintf('%s\\data\\%s\\acc.mat', root, CONT{nt});
        load(fname);
        
        %%
        %AFR
        r1(1)=105; %38 deg
        r2(1)=250; %-35 deg
        c1(1)=323; %-19 deg
        c2(1)=470; %55 deg
        %ASIA
        r1(2)=59; %61 deg
        r2(2)=204; %-12 deg
        c1(2)=475; %57 deg
        c2(2)=720; %180 deg
        %AUS
        r1(3)=201; %-10 deg
        r2(3)=292; %-56 deg
        c1(3)=585; %112 deg
        c2(3)=720; %180 deg
        %CAM
        r1(4)=103; %39 deg
        r2(4)=170; %5 deg
        c1(4)=123; %-119 deg
        c2(4)=240; %-60 deg
        %EUR
        r1(5)=57; %62 deg
        r2(5)=156; %12 deg
        c1(5)=333; %-14 deg
        c2(5)=500; %70 deg
        % NAM
        r1(6)=61; %60 deg
        r2(6)=130; %25 deg
        c1(6)=85; %-138 deg
        c2(6)=256; %-52 deg
        %SAM
        r1(7)=151; %15 deg
        r2(7)=292; %-56 deg
        c1(7)=175; %-93 deg
        c2(7)=296; %-32 deg
        
        for m =1:13
            Rm{m} = R_data{m}(r1(nt):r2(nt),c1(nt):c2(nt));
        end
        
        %     figure(1); clf();imagesc(log(Rm{13})); axis image;
        %     acc=single(acc);
        %     figure(2); clf();imagesc(log(acc)); axis image;colormap(gray);
        
        %%
        % lo-res dimensions:
        
        for m=1:13
            fprintf('m %d ', m);
            
            [nrlo, nclo] = size(Rm{m});
            Nlo = nrlo*nclo;
            
            % hi-res dimensions:
            [nrhi, nchi] = size(acc);
            Nhi = nrhi * nchi;
            
            % size of blocks
            br = nrhi / nrlo;
            bc = nchi / nclo;
            
            % Scaling factor based on ratio of cell areas
            scaling = Nhi / Nlo;
            
            % Create hi-res runoff map
            R = zeros(size(acc),'single');
            
            % fill it with blocks of lo-res runoff
            for r=1:nrlo
                r1 = (r-1)*br+1;
                r2 = r*br;
                for c=1:nclo
                    c1 = (c-1)*bc+1;
                    c2 = c*bc;
                    R(r1:r2, c1:c2) = Rm{m}(r,c) / scaling;
                end;
            end;
            
            %disp('saving')
            matpath = fullfile(root, sprintf('data\\%s\\ISIMIP2B\\R', CONT{nt}));
            if ~isdir(matpath)
                mkdir(matpath);
            end
            
            matfile = fullfile(root, sprintf('data\\%s\\ISIMIP2B\\R\\%s_Rm_%ds_m%d.mat',CONT{nt},fnames(y).name(1:(end-4)),res,m));
            save(matfile,'R','-v7.3');
            
        end
        
        % %%
        % figure(1);
        % imagesc(log(Rm{13}));
        % axis image;
        % colorbar;
        % %%
        %
        % figure(2); clf;
        % imagesc(log(R{13}));
        % axis image;
        % colorbar;
        %
        %
        % %%
        % Rfile = sprintf('Rmonthly_%ds.mat', res);
        % save(Rfile,'R','-v7.3');
    end
end
% end