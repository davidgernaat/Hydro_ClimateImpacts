function ISIMIP_GetNewQs_ISIMIP2B(root,nR,nRCP,sy)
%% Collecting new Qs from ISIMIP downscaled and rerouted discharge maps

if isdeployed
    root2 = 'R:\model\Hydrus';
    nR = str2num(nR);
    nRCP = str2num(nRCP);
    sy = str2num(sy); %Analysis year (2050 or 2100)
end

if ~isdeployed
    clear all
    root = 'W:\share\Hydrus_DOWNQ';
    root2 = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
    root3 = 'R:\model\Hydrus';
    nRCP=2;     %Which RCP
    nR=3;       %Which R-file
end

Continents = {'AFR';'ASIA';'AUS';'CAM';'EUR';'NAM';'SAM'};
scen={'Full','Remain','Eco','Navi'};

RCP_name{1} = 'hist';
RCP_name{2} = 'RCP26';
RCP_name{3} = 'RCP45';
RCP_name{4} = 'RCP60';
RCP_name{5} = 'RCP85';

S=2;    % Which scenario
nC=1;   % Which continent
% m=1;    % Which month
% nR=1;   % Which R-file

%% Load data
disp('Load data')
fname = sprintf('%s\\grid\\Analysis\\Continent_hydropower_Remain_r6_3.mat', root2);
load(fname);

%% Extract file list
disp('Load fnames')
if ~isdeployed
    
    root_fnames = sprintf('%s\\data\\data_prep\\ISIMIP_Runoff\\DB_Runoff_m3s_30y_%s\\%s',root2,RCP_name{nRCP});
    
    fnames_org = dir(root_fnames);
    
    for i=1:numel(fnames_org)
        fnames{i} = fnames_org(i).name(1:end-8);
    end
    
end

%%
if isdeployed
    fname = sprintf('%s\\data\\data_prep\\ISIMIP_fnames\\ISIMIP_fnames_%s.mat', root2, RCP_name{nRCP});
    load(fname);
end

%% Extract new Qs

fprintf('File #%d\n',nR)

for nC=1:7
    fprintf('Continent %s\n',Continents{nC})
    for m=1:13
        fprintf('Month #%d\n',m)
        clear Q
        
        if sy==2100 && nRCP~=1
            fname = sprintf('%s\\%s\\Q_%s_%d_%s_Rm_15s_m%d.mat',root,Continents{nC},RCP_name{nRCP},sy,fnames{nR},m);
        elseif nRCP==1
            fname = sprintf('%s\\%s\\Q_%s_%s_m3s_Rm_15s_m%d.mat',root,Continents{nC},RCP_name{nRCP},fnames{nR},m);
        else
            fname = sprintf('%s\\%s\\Q_%s_%d_%s_m3s_2050_Rm_15s_m%d.mat',root,Continents{nC},RCP_name{nRCP},sy,fnames{nR},m);
        end
        load(fname)
        
        [nr, nc] = size(Q);
        
        continent_in= Continents{nC};
        georef
        
        %% Convert lat/lon to row column and extract Qs
        
        %RD
        for i=1:2000
            if isempty(Basin_COETotRD{nC}{i})==1; QRD{nC}{m}{i}=[]; continue; end;
            clear rRD cRD
            for j=1:numel(Basin_RDlon{nC}{i})
                if isnan(Basin_COETotRD{nC}{i}(j))==1; QRD{nC}{m}{i}(j)=NaN; continue; end;
                
                [rRD(j),cRD(j)] = setpostn(Q,R,Basin_RDlat{nC}{i}(j),Basin_RDlon{nC}{i}(j));
                QRD{nC}{m}{i}(j)= Q(rRD(j),cRD(j));
            end
        end
        
        %P
        for i=1:2000
            if isempty(Basin_COETotP{nC}{i})==1; QP{nC}{m}{i}=[]; continue; end;
            clear rP cP
            for j=1:numel(Basin_Plon{nC}{i})
                if isnan(Basin_COETotP{nC}{i}(j))==1; QP{nC}{m}{i}(j)=NaN; continue; end;
                
                [rP(j),cP(j)] = setpostn(Q,R,Basin_Plat{nC}{i}(j),Basin_Plon{nC}{i}(j));
                QP{nC}{m}{i}(j)= Q(rP(j),cP(j));
            end
        end
        
        
    end
end


%% RD monthly Qs in one vector per dam site
disp('Put months in one vector RD')
for nC=1:7
    for i=1:2000
        
        if isempty(QRD{nC}{m}{i})==1; QRDm{nC}{i}=[]; continue; end;
        
        for j=1:numel(QRD{nC}{1}{i})
            for m=1:13
                QRDm{nC}{i}{j}(m) = QRD{nC}{m}{i}(j);
            end
        end
    end
end

%% P monthly Qs in one vector per dam site
disp('Put months in one vector P')
for nC=1:7
    for i=1:2000
        
        if isempty(QP{nC}{m}{i})==1; QPm{nC}{i}=[]; continue; end;
        
        for j=1:numel(QP{nC}{1}{i})
            for m=1:13
                QPm{nC}{i}{j}(m) = QP{nC}{m}{i}(j);
            end
        end
    end
end

%% Calculate Q_design, Q_design_LF and Q_design_mean. See Monthly_loadfactor_output.m
disp('Calculate Q_designs RD')
%RD
for nC=1:7
    for i=1:2000
        
        if isempty(QRDm{nC}{i})==1; QdesignRD{nC}{i} = []; LF_QdesignRD{nC}{i} = [];
            Qdesign_meanRD{nC}{i} = []; QmeanRD{nC}{i} = []; continue; end;
        
        for j=1:numel(QRDm{nC}{i})
            clear LF_Qdesign_Monthly
            if isnan(QRDm{nC}{i}{j}(1))==1; QdesignRD{nC}{i}(j) = NaN; LF_QdesignRD{nC}{i}(j) = NaN;
                Qdesign_meanRD{nC}{i}(j) = NaN; QmeanRD{nC}{i}(j) = NaN; continue; end;
            
            %Q_design
            Qsort{nC}{i}{j} = sort(QRDm{nC}{i}{j},'descend');
            QdesignRD{nC}{i}(j) = Qsort{nC}{i}{j}(4); %forth highest discharge month
            
            %Q_design_LF
            for m=1:12; LF_Qdesign_Monthly(m) = min(1,Qsort{nC}{i}{j}(m)/QdesignRD{nC}{i}(j))'; end
            LF_QdesignRD{nC}{i}(j) = mean(LF_Qdesign_Monthly);
            
            %Q_design_mean
            Qdesign_meanRD{nC}{i}(j) = mean(min(QdesignRD{nC}{i}(j),Qsort{nC}{i}{j})); %Average flow rate given Qdesign
            
            %Qoutlets (for seasonal storage systems)
            QmeanRD{nC}{i}(j) = QRDm{nC}{i}{j}(13);
            
        end
    end
end

%     nL=2366;
%     nL=494;
%     nBB=10;
%     fprintf('QdesignRD ISIMIP = %0.2f\n',QdesignRD{nC}{nBB}(nL))
%     fprintf('QdesignRD LPJmL = %0.2f\n\n',Basin_Qoutlets_design{nC}{nBB}(nL))

%%
disp('Calculate Q_designs P')
clear Qsort
%P
for nC=1:7
    for i=1:2000
        
        if isempty(QPm{nC}{i})==1; QdesignP{nC}{i} = []; LF_QdesignP{nC}{i} = [];
            Qdesign_meanP{nC}{i} = []; QmeanP{nC}{i} = []; continue; end;
        
        for j=1:numel(QPm{nC}{i})
            clear LF_Qdesign_Monthly
            if isnan(QPm{nC}{i}{j}(1))==1; QdesignP{nC}{i}(j) = NaN; LF_QdesignP{nC}{i}(j) = NaN;
                Qdesign_meanP{nC}{i}(j) = NaN; QmeanP{nC}{i}(j) = NaN; continue; end;
            
            %Q_design
            Qsort{nC}{i}{j} = sort(QPm{nC}{i}{j},'descend');
            QdesignP{nC}{i}(j) = Qsort{nC}{i}{j}(4); %forth highest discharge month
            
            %Q_design_LF
            for m=1:12; LF_Qdesign_Monthly(m) = min(1,Qsort{nC}{i}{j}(m)/QdesignP{nC}{i}(j))'; end
            LF_QdesignP{nC}{i}(j) = mean(LF_Qdesign_Monthly);
            
            %Q_design_mean
            Qdesign_meanP{nC}{i}(j) = mean(min(QdesignP{nC}{i}(j),Qsort{nC}{i}{j})); %Average flow rate given Qdesign
            
            %Qoutlets (for seasonal storage systems)
            QmeanP{nC}{i}(j) = QPm{nC}{i}{j}(13);
        end
    end
end


%% Save
if nRCP~=1
    matfile = fullfile(root, sprintf('ISIMIP_NewQs_Remain_%s_%d_%s.mat',RCP_name{nRCP},sy,fnames{nR}));
else
    matfile = fullfile(root, sprintf('ISIMIP_NewQs_Remain_%s_%s.mat',RCP_name{nRCP},fnames{nR}));
end

save(matfile,'-v7.3','QdesignRD','LF_QdesignRD','Qdesign_meanRD','QmeanRD','QdesignP','LF_QdesignP','Qdesign_meanP','QmeanP');

