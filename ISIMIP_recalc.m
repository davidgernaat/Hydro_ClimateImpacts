%% ISIMIP re-calculate energy potential and costs
% Step 1: re-collect Qs based on GCM and GHM runoff maps from ISIMIP. See ISIMIP_GetNewQs.m
% Step 2: re-run cost models that calculate energy potential and costs.

clear all
root2 = 'Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\Global';
root = 'W:\share\Hydrus_DOWNQ';
Continents = {'AFR';'ASIA';'AUS';'CAM';'EUR';'NAM';'SAM'};
scen={'Full','Remain','Eco','Navi'};
addpathname = sprintf('%s\\Hydrus\\functions', root2);
addpath(addpathname)
output=1; %1 = output / 0 = no output
nc=[1 2 3 4 5 6 7];

RCP_name{1} = 'hist';
RCP_name{2} = 'RCP26';
RCP_name{3} = 'RCP45';
RCP_name{4} = 'RCP60';
RCP_name{5} = 'RCP85';
RCP_name{6} = 'ISIMIP2B';

pv = {'Exis','Remain'};

S=2;    % Which scenario
nC=1;   % Which continent
m=1;    % Which month
nR=3;   % Which R-file
ir=6;   % Which RCP level
sy = 2100; % Which study year
p = 2;  % Which potential (existing or remaining)

%% Load data
disp('Load data')
fname = sprintf('%s\\grid\\Analysis\\Continent_hydropower_%s_r6_5.mat', root2, scen{S});
load(fname);

%% The ISMIP file names
% for ir=2:5
    clear R_names root_data fnames
    
    if ir==1
        root_data = sprintf('%s\\data\\data_prep\\ISIMIP_NewQs\\%s\\%s',root2,RCP_name{ir},pv{p});
    elseif ir==6
        root_data = sprintf('%s\\data\\data_prep\\ISIMIP_NewQs_2B',root2);
    else
        root_data = sprintf('%s\\data\\data_prep\\ISIMIP_NewQs\\%s_%d\\%s',root2,RCP_name{ir},sy,pv{p});
    end
    
    fnames = dir(root_data);
    
    for j=1:numel(fnames)
        R_names{j} = sprintf('%s\\%s',root_data,(fnames(j).name));
    end
    
    %% Filling output variables
    % for C=1:7
    %     for i=1:2000
    %         for j=1:numel(Basin_COETotRD{i}{j})
    %             RDPnet{C}{i}(j)=0;
    %             COETotRD{C}{i}(j)=0;
    %         end
    %     end
    % end
    % for C=1:7
    %     for i=1:2000
    %         for j=1:4
    %             if isempty(Basin_COETotP{C}{i})==1
    %                 PPnet{C}{i}{j}=[];
    %                 COEP{C}{i}{j}=[];
    %             else
    %                 for k=1:numel(Basin_COETotP{C}{i}{j})
    %                     PPnet{C}{i}{j}(l)=0;
    %                     COEP{C}{i}{j}(l)=0;
    %                 end
    %             end
    %         end
    %     end
    % end
    %% Re calculation
    for nR=3:numel(R_names)
        clear Basin_COETotPs Basin_COETotRDs Basin_PPnets Basin_RDPnets QdesignPA Qdesign_meanPA LF_QdesignPA
        clear Basin_PPnetsA Basin_COETotPsA COETotRD RDPnet PPnet COEP
        
        % nR=11;
        fprintf('Loading %s NewQs: %d of %d\n',RCP_name{ir} ,nR, numel(R_names))
        load(R_names{nR});
        
        % fprintf('Loading direct Hydrus output\n')
        % fname = sprintf('%s\\output\\Remain\\NAM\\Basin4_output.mat', root2);
        % load(fname);
        
        %% Vector extending
        for C=1:7
            for B=1:2000
                for k=(numel(Basin_RDDepth{C}{B})+1):numel(Basin_COETotRD{C}{B})
                    Basin_RDDepth{C}{B}(k)=0;
                    Basin_DisOutlet{C}{B}(k)=0;
                    Basin_Qoutlets_design_LF{C}{B}(k)=0;
                    Basin_Qoutlets_design{C}{B}(k)=0;
                    Basin_Zoutlets{C}{B}(k)=0;
                end
                for k=(numel(QdesignRD{C}{B})+1):numel(Basin_COETotRD{C}{B})
                    QdesignRD{C}{B}(k)=0;
                    LF_QdesignRD{C}{B}(k)=0;
                end
            end
        end
        
        %% Find indices of seperate dam systems in PnetAlls
        for C=1:7
            for B=1:2000
                PID{C}{B}  = find(Basin_SysID{C}{B}==1);
                RDID{C}{B} = find(Basin_SysID{C}{B}==2);
            end
        end
        
        %% Create vars with deselected dam info
        for C=1:7
            for B=1:2000
                Basin_COETotPs{C}{B} = Basin_COEAlls{C}{B}(PID{C}{B});
                Basin_COETotRDs{C}{B} = Basin_COEAlls{C}{B}(RDID{C}{B});
                
                Basin_PPnets{C}{B} = Basin_PnetAlls{C}{B}(PID{C}{B});
                Basin_RDPnets{C}{B} = Basin_PnetAlls{C}{B}(RDID{C}{B});
            end
        end
        
        %%  RD costs model
        disp('RD costmodel')
        for C=1:7
            fprintf('Continent #%d of #7\n',C)
            for B=1:2000
                
                if isempty(Basin_COETotRDs{C}{B})==1;
                    COETotRD{C}{B}=[]; RDPnet{C}{B}=[]; COETotRD_HYD{C}{B}=[]; RDPnet_HYD{C}{B}=[]; OptSpecCapRD{C}{B}=[];
                    continue; end
                
                for i=1:numel(Basin_COETotRDs{C}{B})
                    
                    if isnan(Basin_COETotRDs{C}{B}(i))==1
                        COETotRD_HYD{C}{B}(i)=NaN; RDPnet_HYD{C}{B}(i)=NaN;
                        COETotRD{C}{B}(i)=NaN; RDPnet{C}{B}(i)=NaN; OptSpecCapRD{C}{B}(i)=NaN; continue; end;
                    
                    % As in Hydrus-model
                    [COETotRD_HYD{C}{B}(i), RDPnet_HYD{C}{B}(i), ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                        RDcostmodel(Basin_OptDL{C}{B}(i),Basin_OptDH{C}{B}(i),Basin_OptPop{C}{B}(i),Basin_Qoutlets_design{C}{B}(i), ...
                        Basin_Qoutlets_design_LF{C}{B}(i),Basin_DisOutlet{C}{B}(i),Basin_OptLV{C}{B}(i),Basin_RDDepth{C}{B}(i),1,1,1,10,0);
                    
                    RDPnet_HYD{C}{B}(i) = RDPnet_HYD{C}{B}(i)*1e-6; % GWh
                    
                    % As with ISIMIP NewQs
                    [COETotRD{C}{B}(i), RDPnet{C}{B}(i),~, ~, ~, ~, ~, ~, ~, ~, OptSpecCapRD{C}{B}(i)] = ...
                        RDcostmodel(Basin_OptDL{C}{B}(i),Basin_OptDH{C}{B}(i),Basin_OptPop{C}{B}(i),...
                        QdesignRD{C}{B}(i), LF_QdesignRD{C}{B}(i),... %New data
                        Basin_DisOutlet{C}{B}(i),Basin_OptLV{C}{B}(i),Basin_RDDepth{C}{B}(i),1,1,0,10,0);
                    
                    RDPnet{C}{B}(i) = RDPnet{C}{B}(i)*1e-6; % GWh
                    
                    
                end
            end
        end
        
        C=2;
        B=1;
                fprintf('Direct Hydrus RD potential: %0.0f GWh\n',sum(Basin_RDPnets{C}{B}(~isnan(Basin_RDPnets{C}{B}(:)))));
                fprintf('Re-calculated RD potential: %0.0f GWh\n',sum(RDPnet_HYD{C}{B}(~isnan(RDPnet_HYD{C}{B}(:)))));
%                 fprintf('Re-calculated RD potential corrected: %0.0f GWh\n',sum(RDPnet_HYD{11}{C}{B}(~isnan(Basin_RDPnets{C}{B}(:)))));
        %
                fprintf('ISIMIP RD potential: %0.0f GWh\n',sum(RDPnet{C}{B}(~isnan(RDPnet{C}{B}(:)))));
%                 fprintf('ISIMIP RD potential: %0.0f GWh\n',sum(RDPnet{11}{C}{B}(~isnan(Basin_RDPnets{C}{B}(:)))));
        
        %% Vector extending
        for C=1:7
            for B=1:2000
                for k=(numel(QdesignP{C}{B})+1):numel(Basin_COETotP{C}{B})
                    QdesignP{C}{B}(k)=0;
                    Qdesign_meanP{C}{B}(k)=0;
                    LF_QdesignP{C}{B}(k)=0;
                end
            end
        end
        
        %% Make array of P-variables to get them in the required loops
        disp('Prepping data arrays for P cost model')
        
        for C=1:7
            for B=1:2000
                n=numel(Basin_COETotRD{C}{B});
                a = [n+1;(2*n)+1;(3*n+1)];
                j=1;
                k=0;
                for i=1:numel(Basin_COETotP{C}{B})
                    if ismember(i,a); j=j+1; k=0; end
                    k=k+1;
                    
                    QdesignPA{C}{B}{j}(k) = QdesignP{C}{B}(i);
                    Qdesign_meanPA{C}{B}{j}(k) = Qdesign_meanP{C}{B}(i);
                    LF_QdesignPA{C}{B}{j}(k) = LF_QdesignP{C}{B}(i);
                    
                    Basin_PPnetsA{C}{B}{j}(k) = Basin_PPnets{C}{B}(i);
                    Basin_COETotPsA{C}{B}{j}(k) = Basin_COETotPs{C}{B}(i);
                    
                end
            end
        end
        
        %% P Cost
        disp('P cost model')
        for C=1:7
            fprintf('Continent #%d of #7\n',C)
            for B=1:2000
                
                if isempty(Basin_ZPinlet{C}{B})==1;
                    COEP{C}{B}=[]; PPnet{C}{B}=[]; COEP_HYD{C}{B}=[]; PPnet_HYD{C}{B}=[]; OptSpecCapP{C}{B}=[];
                    continue; end
                
                for j=1:4
                    for i=1:numel(Basin_ZPinlet{C}{B}{j})
                        
                        if isnan(Basin_COETotPsA{C}{B}{j}(i))==1;
                            COEP{C}{B}{j}(i)=NaN; PPnet{C}{B}{j}(i)=NaN; COEP_HYD{C}{B}{j}(i)=NaN; PPnet_HYD{C}{B}{j}(i)=NaN;
                            OptSpecCapP{C}{B}{j}(i)=NaN;
                            continue; end;
                        
                        % As in Hydrus-model
                        [COEP_HYD{C}{B}{j}(i), PPnet_HYD{C}{B}{j}(i), ~, ~, ~, ~, ~] = ...
                            costmodel_pipesys(Basin_Zoutlets{C}{B}(i),Basin_ZPinlet{C}{B}{j}(i),Basin_PL{C}{B}{j}(i),Basin_QDesignPinlet{C}{B}{j}(i), ...
                            Basin_QDesignMeanPinlet{C}{B}{j}(i), Basin_QDesignLFPinlet{C}{B}{j}(i),Basin_DisOutlet{C}{B}(i),0,1000,0);
                        
                        PPnet_HYD{C}{B}{j}(i) = PPnet_HYD{C}{B}{j}(i)*1e-6; %GWh
                        
                        % As with ISIMIP NewQs
                        [COEP{C}{B}{j}(i), PPnet{C}{B}{j}(i), ~, ~, ~, ~, OptSpecCapP{C}{B}{j}(i), ~] = ...
                            costmodel_pipesys(Basin_Zoutlets{C}{B}(i),Basin_ZPinlet{C}{B}{j}(i),Basin_PL{C}{B}{j}(i),...
                            QdesignPA{C}{B}{j}(i), Qdesign_meanPA{C}{B}{j}(i), LF_QdesignPA{C}{B}{j}(i),...
                            Basin_DisOutlet{C}{B}(i),0,1000,0);
                        
                        PPnet{C}{B}{j}(i) = PPnet{C}{B}{j}(i)*1e-6; %GWh
                        
                    end
                end
            end
        end
        
        C=1;
        B=1;
        fprintf('Direct Hydrus P potential: %0.0f GWh\n',sum(Basin_PPnets{C}{B}(~isnan(Basin_PPnets{C}{B}(:)))));
        PPnet_HYDa{C}{B}  = horzcat(PPnet_HYD{C}{B}{:})';
        fprintf('Re-calculated P potential: %0.0f GWh\n',sum(PPnet_HYDa{C}{B}(~isnan(PPnet_HYDa{C}{B}(:)))));
        % fprintf('Re-calculated P potential corrected: %0.0f GWh\n',sum(PPnet_HYDa{11}{C}{B}(~isnan(Basin_PPnets{C}{B}(:)))));
        %
        PPneta{C}{B}  = horzcat(PPnet{C}{B}{:})';
        fprintf('ISIMIP P potential: %0.0f GWh\n',sum(PPneta{C}{B}(~isnan(PPneta{C}{B}(:)))));
        % fprintf('ISIMIP P potential corrected: %0.0f GWh\n\n',sum(PPneta{11}{C}{B}(~isnan(Basin_PPnets{C}{B}(:)))));
        
        
        %%
        %     %% Check RD
        %     nR=11;
        %     C=1;
        %     B=10;
        %     i=494;
        %     %         i=490;
        %     %         i=1211;
        %     fprintf('RD\n');
        %     fprintf('Hydrus Pnet = %0.2f\n',Basin_RDPnet{C}{B}(i));
        %     fprintf('ISIMIP Pnet = %0.2f\n\n',RDPnet{nR}{C}{B}(i));
        %     fprintf(' Hydrus COE = %0.4f\n',Basin_COETotRD{C}{B}(i));
        %     fprintf(' ISIMIP COE = %0.4f\n\n',COETotRD{nR}{C}{B}(i));
        %
        %     %% Check P
        %     nR=11;
        %     C=1;
        %     B=10;
        %     j=1;
        %     i=276;
        %     fprintf('P\n')
        %     fprintf('Hydrus Pnet = %0.2f\n',Basin_PPneta{C}{B}{j}(i));
        %     fprintf('ISIMIP Pnet = %0.2f\n\n',PPnet{nR}{C}{B}{j}(i));
        %     fprintf(' Hydrus COE = %0.4f\n',Basin_COETotPa{C}{B}{j}(i));
        %     fprintf(' ISIMIP COE = %0.4f\n',COEP{nR}{C}{B}{j}(i));
        
        % Save
        if output==1
            disp('Save')
            
            if ir==1
                matpath = fullfile(root2, sprintf('output\\ISIMIP\\%s\\%s',RCP_name{ir},pv{p}));
            elseif ir==6
                matpath = fullfile(root2, sprintf('output\\%s',RCP_name{ir}));
            else
                matpath = fullfile(root2, sprintf('output\\ISIMIP\\%s_%d\\%s',RCP_name{ir},sy,pv{p}));
            end
            
            if ~isdir(matpath)
                mkdir(matpath);
            end
            
            if ir==1
                matfile = fullfile(root2, sprintf('output\\ISIMIP\\%s\\%s\\ISIMIP_NewPnet_%s_%s_%s',RCP_name{ir},pv{p},RCP_name{ir},pv{p},fnames(nR).name));
            elseif ir==6
                matfile = fullfile(root2, sprintf('output\\%s\\ISIMIP2B_NewPnet_%s',RCP_name{ir},fnames(nR).name(35:end)));
            else
                matfile = fullfile(root2, sprintf('output\\ISIMIP\\%s_%d\\%s\\ISIMIP_NewPnet_%s_%d_%s_%s',RCP_name{ir},sy,pv{p},RCP_name{ir},sy,pv{p},fnames(nR).name));
            end
            save(matfile,'-v7.3','RDPnet','COETotRD','OptSpecCapRD','PPnet','COEP','OptSpecCapP');
        end
        
    end
%     end

%% Check RD potentials
% % nR=11;
% % C=7;
% % B=1;
% % Continents = {'AFR';'ASIA';'AUS';'CAM';'EUR';'NAM';'SAM'};
% %
% for C=1:7
%     for B=1:10 %numel(RDPnet{nR}{C})
%         Tot_HYD{C}(B) = sum(RDPnet_HYD{nR}{C}{B}(find(~isnan(RDPnet_HYD{nR}{C}{B}(:)))));
%         Tot_ISI{C}(B) = sum(RDPnet{nR}{C}{B}(find(~isnan(RDPnet{nR}{C}{B}(:)))));
%
%         %     fprintf('Basin #%d\n',B);
%         %     fprintf('Hydrus Pnet = %0.2f TWh\n',Tot_HYD{C}(B) *1e-3);
%         %     fprintf('ISIMIP Pnet = %0.2f TWh\n\n',Tot_ISI{C}(B) *1e-3);
%     end
% %     fprintf('%s Hydrus Pnet = %0.2f TWh\n',Continents{C}, sum(Tot_HYD{C}(:)) *1e-3);
% %     fprintf('%s ISIMIP Pnet = %0.2f TWh\n\n',Continents{C}, sum(Tot_ISI{C}(:)) *1e-3);
%
% end
%
% % fprintf('Total Hydrus Pnet = %0.2f TWh\n',sum(Tot_HYD(:)) *1e-3);
% % fprintf('Total ISIMIP Pnet = %0.2f TWh\n\n',sum(Tot_ISI(:)) *1e-3);
%
% Tot_HYDW2 = horzcat(Tot_HYD{:});
% Tot_ISIW2 = horzcat(Tot_ISI{:});
% TotHYDW = sum(Tot_HYDW2(:));
% TotISIW = sum(Tot_ISIW2(:));
%
% fprintf('World Hydrus RDPnet = %0.2f TWh\n',TotHYDW *1e-3);
% fprintf('World ISIMIP RDPnet = %0.2f TWh\n\n',TotISIW *1e-3);
%
%
% %% Check P potentials
% % nR=11;
% C=7;
% B=1;
% Continents = {'AFR';'ASIA';'AUS';'CAM';'EUR';'NAM';'SAM'};
%
% for C=1:7
%     for B=1:10
%
%         Tot_HYDP2{C}{B} = horzcat(PPnet_HYD{nR}{C}{B}{:});
%         Tot_ISIP2{C}{B} = horzcat(PPnet{nR}{C}{B}{:});
%
%         Tot_HYD{C}(B) = sum(Tot_HYDP2{C}{B}(find(~isnan(Tot_HYDP2{C}{B}(:)))));
%         Tot_ISI{C}(B) = sum(Tot_ISIP2{C}{B}(find(~isnan(Tot_ISIP2{C}{B}(:)))));
%
%         %     fprintf('Basin #%d\n',B);
%         %     fprintf('Hydrus Pnet = %0.2f TWh\n',Tot_HYD{C}(B) *1e-3);
%         %     fprintf('ISIMIP Pnet = %0.2f TWh\n\n',Tot_ISI{C}(B) *1e-3);
%     end
% %     fprintf('%s Hydrus Pnet = %0.2f TWh\n',Continents{C}, sum(Tot_HYD{C}(:)) *1e-3);
% %     fprintf('%s ISIMIP Pnet = %0.2f TWh\n\n',Continents{C}, sum(Tot_ISI{C}(:)) *1e-3);
%
% end
%
% % fprintf('Total Hydrus Pnet = %0.2f TWh\n',sum(Tot_HYD(:)) *1e-3);
% % fprintf('Total ISIMIP Pnet = %0.2f TWh\n\n',sum(Tot_ISI(:)) *1e-3);
%
% Tot_HYDW2 = horzcat(Tot_HYD{:});
% Tot_ISIW2 = horzcat(Tot_ISI{:});
% TotHYDW = sum(Tot_HYDW2(:));
% TotISIW = sum(Tot_ISIW2(:));
%
% fprintf('World Hydrus  PPnet = %0.2f TWh\n',TotHYDW *1e-3);
% fprintf('World ISIMIP  PPnet = %0.2f TWh\n\n',TotISIW *1e-3);
%
% %%
% disp('Load data')
% fname = sprintf('%s\\grid\\Analysis\\ISIMIP_NewPnet_%s.mat', root2, scen{S});
% load(fname);
%
%
% %% Check if Hydrus output matches the output aggregated
% %Result is eactly the same


%% load Basin4_output.mat
% Direct Hydrus output
% aDPPnetEnd4  = horzcat(PPnetend{:})'; % GWh  Dam-Pipe systems Pnet
% fprintf('\nTot PnetAlls potential: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
% fprintf('Tot P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
% fprintf('Tot RD potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet))));
% fprintf('Tot RD + P potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet)))+sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
%
% %Aggregared Hydrus output via grid output
% fprintf('\nTot PnetAlls potential: %0.0f GWh\n',sum(Basin_PnetAlls{6}{4}(~isnan(Basin_PnetAlls{6}{4}(:)))));
% fprintf('Tot P potential: %0.0f GWh\n',sum(Basin_PPnet{6}{4}(~isnan(Basin_PPnet{6}{4}(:)))));
% fprintf('Tot RD potential: %0.0f GWh\n',sum(Basin_RDPnet{6}{4}(~isnan(Basin_RDPnet{6}{4}(:)))));
% fprintf('Tot RD + P potential: %0.0f GWh\n',sum(Basin_PPnet{6}{4}(~isnan(Basin_PPnet{6}{4}(:))))....
%     +sum(Basin_RDPnet{6}{4}(~isnan(Basin_RDPnet{6}{4}(:)))));

% % Re-calculated Hydrus output
% % fprintf('\nTot PnetAlls potential: %0.0f GWh\n',sum(PnetAlls(~isnan(PnetAlls))));
% % fprintf('Tot P potential: %0.0f GWh\n',sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));
% fprintf('Tot RD potential: %0.0f GWh\n',sum(RDPnet_HYD{11}{6}{4}(~isnan(RDPnet_HYD{11}{6}{4}(:)))));
% % fprintf('Tot RD + P potential: %0.0f GWh\n',sum(RDPnet(~isnan(RDPnet)))+sum(aDPPnetEnd4(~isnan(aDPPnetEnd4))));

%% Check Global Potential direct
% 
% for C=1:7
%     RDPnet_HYD_W1{C} = horzcat(RDPnet_HYD{C}{:});
%     for B=1:2000
%         if isempty(PPnet_HYD{C}{B}); PPnet_HYD_W2{C}{B}=[]; continue; end;
%         PPnet_HYD_W2{C}{B} = horzcat(PPnet_HYD{C}{B}{:});
%     end
%     PPnet_HYD_W1{C} = horzcat(PPnet_HYD_W2{C}{:});
% end
% 
% RDPnet_HYD_W = horzcat(RDPnet_HYD_W1{:});
% PPnet_HYD_W = horzcat(PPnet_HYD_W1{:});
% 
% RDPnet_HYD_Wsum = sum(RDPnet_HYD_W(find(~isnan(RDPnet_HYD_W(:)))));
% PPnet_HYD_Wsum = sum(PPnet_HYD_W(find(~isnan(PPnet_HYD_W(:)))));
% 
% Pnet_HYD_Wsum = RDPnet_HYD_Wsum*1e-6 + PPnet_HYD_Wsum;
% 
% fprintf('Global Potential = %0.2f PWh\n',Pnet_HYD_Wsum*1e-6);
% fprintf('Global Potential RD = %0.2f PWh\n',RDPnet_HYD_Wsum*1e-6*1e-6);
% fprintf('Global Potential P  = %0.2f PWh\n\n',PPnet_HYD_Wsum*1e-6);

%% Check Global Potential recalculated based on direct variables
% 
% for C=1:7
%     RDPnet_W1{11}{C} = vertcat(Basin_RDPnets{C}{:});
%     PPnet_W2{11}{C} = vertcat(Basin_PPnets{C}{:});
% end
% 
% RDPnet_W{11} = vertcat(RDPnet_W1{11}{:});
% PPnet_W{11} = horzcat(PPnet_W1{11}{:});
% 
% RDPnet_Wsum{11} = sum(RDPnet_W{11}(find(~isnan(RDPnet_W{11}(:)))));
% PPnet_Wsum{11} = sum(PPnet_W{11}(find(~isnan(PPnet_W{11}(:)))));
% 
% Pnet_Wsum{11} = RDPnet_Wsum{11} + PPnet_Wsum{11};
% 
% fprintf('Global Potential = %0.2f PWh\n',Pnet_Wsum{11}*1e-6);
% fprintf('Global Potential RD = %0.2f PWh\n',RDPnet_Wsum{11}*1e-6);
% fprintf('Global Potential P  = %0.2f PWh\n\n',PPnet_Wsum{11}*1e-6);

%% Check Global Potential recalculated based on ISIMIP
% 
% for C=1:7
%     RDPnet_HYD_W1{C} = horzcat(RDPnet{C}{:});
%     for B=1:2000
%         if isempty(PPnet{C}{B}); PPnet_HYD_W2{C}{B}=[]; continue; end;
%         PPnet_HYD_W2{C}{B} = horzcat(PPnet{C}{B}{:});
%     end
%     PPnet_HYD_W1{C} = horzcat(PPnet_HYD_W2{C}{:});
% end
% 
% RDPnet_HYD_W = horzcat(RDPnet_HYD_W1{:});
% PPnet_HYD_W = horzcat(PPnet_HYD_W1{:});
% 
% RDPnet_HYD_Wsum = sum(RDPnet_HYD_W(find(~isnan(RDPnet_HYD_W(:)))));
% PPnet_HYD_Wsum = sum(PPnet_HYD_W(find(~isnan(PPnet_HYD_W(:)))));
% 
% Pnet_HYD_Wsum = RDPnet_HYD_Wsum + PPnet_HYD_Wsum;
% 
% fprintf('Global Potential = %0.2f PWh\n',Pnet_HYD_Wsum*1e-6);
% fprintf('Global Potential RD = %0.2f PWh\n',RDPnet_HYD_Wsum*1e-6);
% fprintf('Global Potential P  = %0.2f PWh\n\n',PPnet_HYD_Wsum*1e-6);

%% Continental check
% for C=1:7
%     %PPnet_HYD_W1C{11}{C} = sum(PPnet_HYD_W1{11}{C}(find(~isnan(PPnet_HYD_W1{11}{C}(:)))));
%     %fprintf('P Potential %s = %0.02f\n',Continents{C}, PPnet_HYD_W1C{11}{C}*1e-3)
%
%     RDPnet_HYD_W1C{11}{C} = sum(RDPnet_HYD_W1{11}{C}(find(~isnan(RDPnet_HYD_W1{11}{C}(:)))));
%     fprintf('RD Potential %s = %0.02f\n',Continents{C}, RDPnet_HYD_W1C{11}{C}*1e-3)
% end
%
% %%
% C=1;
% RDPnetAFR_Calc = horzcat(RDPnet_HYD{11}{C}{:})';
% RDPnetAFR_Hydrus  = vertcat(Basin_RDPnets{C}{:});
%
% RDPnetAFR_CalcSum = sum(RDPnetAFR_Calc(find(~isnan(RDPnetAFR_Calc(:)))));
% RDPnetAFR_HydrusSum = sum(RDPnetAFR_Hydrus(find(~isnan(RDPnetAFR_Hydrus(:)))));
%
% fprintf('Calc Potential = %0.2f\n',RDPnetAFR_CalcSum);
% fprintf('hydrus Potential = %0.2f\n',RDPnetAFR_HydrusSum);

%% Check for weird numbers RD
% for C=1:7
%         COETotRD2{C} = horzcat(COETotRD{C}{:});
% end
% COETotRD3 = horzcat(COETotRD2{:});
% a=tabulate(COETotRD3(:));

%% Check for weird numbers P
% for C=1:7
%     for B=1:2000
%         if isempty(COEP{C}{B})==1; continue; end
%         COEP2{C}{B} = horzcat(COEP{C}{B}{:});
%     end
%     COEP3{C} = horzcat(COEP2{C}{:});
% end
% COEP4 = horzcat(COEP3{:});
% b=tabulate(COEP4(:));

%% save 
% matfile = fullfile(root2, sprintf('grid\\Analysis\\Continent_hydropower_Remain_r6_5_latlon.mat'));
% save(matfile,'-v7.3','Basin_COEAlls','Basin_PnetAlls','Basin_COETotP','Basin_COETotRD','Basin_Plat','Basin_Plon','Basin_RDlat','Basin_RDlon','RegionID','CountryID','Basin_SysID');