function [COEminRoR PnetRoRmin CostElements npipemin OptInv OptP OptSpecCap Opthf OptDP] = costmodel_pipesys(Zout,Zin,L,Qdesign,Qdesign_mean,LF,Dis,cost_constr,cost_lim,Quakerate)

%% Variables from outer space
% clear all
% figure(1); clf;
% figure(2); clf;
% figure(3); clf;

% k= 5938;
% l=1;
% m=1;
% 
% Zout   = 96;%Zoutlets(k);
% Zin     = 377;%ZPinlet{k}{l}(m);
% L       = 1753.34;%PL{k}{l}(m);       %Pipe lengt (m)
% Qdesign = 243;%QDesignPinlet{k}{l}(m);
% Qdesign_mean = 110.12;%QDesignMeanPinlet{k}{l}(m);
% LF = 0.49;%QDesignLFPinlet{k}{l}(m);
% Dis = 17;%DisOutlet(k);
% cost_constr=0;
% cost_lim=1000;
% Quakerate=0;

% Qmd = zeros(12,1);
% for o=1:12; Qmd(o) = Qm{o}(inlet{j}{m}(k)); end;

% Qdesign is allows 30% exceedence, so 70% of maximum flow
% Here we take the 4th highest mean monthly discharge
% Qsort = sort(Qmd,'descend');
% Qdesign = Qsort(4); %forth highest discharge month
% Qdesign_mean = mean(min(Qdesign,Qsort)); %Average flow rate given Qdesign

% Loadfactor is calculated by deviding montly discharge by Qdesign.
% If higher than LF=1.
% for i=1:12; LF_Qdesign_Monthly(i) = min(1,Qsort(i)/Qdesign)'; end
% LF_Qdesign = mean(LF_Qdesign_Monthly);

%% constants
npipe       = linspace(1,4,4);  % Number of pipes
Hgross      = Zin-Zout;        %Gross head (m)
D           = linspace(0.1,4,10);  %Pipe diameter (meters)
i           = 0;
cc          = hsv(4);         %Color generator
g           = 9.8;             %Gravitational acceleration (m/s2)
rho         = 1000;            %Density of water (kg/m3)
mu          = 0.001;           %Fluid viscosity of water (N-s/m2))
e           = 0.2;             %Roughness constant rough concrete (m)
eta         = 0.8;             %Efficiency turbine
lifetime    = 40;              %Years
syear       = (60*60*24*365);  %Seconds in year
V           = 1.25 * 400e6;    %Volume of lake (M3)
%Qtot        = V/syear;        %Available flow by Dixcence lake
p           = 0.05;            %Price per kWh ($/kWh)
interest    = 0.1;             %Interest on capital
DW          = 50;              %Dam width (m)
OMshare     = 0.0225;          %As share of total investments
ER          = 0.1725;          %Exchange rate NOK to US$ 01/01/2010
IR2002      = 1.24;            %Inflation rate conversion t0 2010 USD: World Bank Real effective exchange rate index (2010 = 100)
IR2005      = 1.09;            %Inflation rate conversion t0 2010 USD: http://data.worldbank.org/indicator/PX.REX.REER?page=2
IR2013      = 0.99;            %Inflation rate conversion to 2010 USD: http://data.worldbank.org/indicator/PX.REX.REER
TransCost   = 0.0034 * IR2013; %Fixed transmission cost per delivered kWh: http://www.eia.gov/forecasts/aeo/tables_ref.cfm Table 8 EIA ANNUAL ENERGY OUTLOOK 2015
Ownersrate  = 0.25;            %Owners cost due to lead times IRENA 2012 Hydropower

%% Physical calculations

for j=1:numel(npipe)
    
    Atunnel = pi*(D/2).^2;   %Cross sectional Tunnel(m2)
    f = 0.01;                %Initial friction number
    Q = Qdesign;
    u = (Qdesign/npipe(j))./Atunnel; % Water velocity (m/s)
    
    % Calculate friction factor in iterations
    for c = 1:10
        hf = f*L.*u.^2 ./ (2*g*D);                            % Head loss
        Re = D.*rho.*u./mu;                                   % Reynold number
        RHS = -2.*log10(e./(3.71.*D) + 2.51./(Re.*sqrt(f)));  % RHS of Colebrook-White eq
        f = (1./RHS).^2;                                      % New friction factor
    end
    
    for i=1:numel(D)
        if hf(i)>Hgross; hf(i)=NaN; end
    end
    
    PTgross = rho*g*Qdesign_mean*Hgross*8760*1e-9;    % Theoretical gross potential
    PTnet = rho*g*Qdesign_mean*(Hgross-hf)*8760*1e-9; % Theoretical net potential
    Pror = rho*g*(Hgross-hf).*eta*Qdesign;            % Capacity based on calculated Q (W)
    Pror_end{j} = Pror;
    hf_end{j} = hf/Hgross;                            % Percentage loss due to head loss
    PnetRoR = (Pror*1e-3*8760).*LF;               % Yearly production with LF (kWh)
    
    %% Cost calculations
    
    %Different cost metrics based on capacity, Q, head and dam height
    %CostGunduz = -0.0254 + 1.81*QQ + 0.0848*(Hgross-Ahf) + 0.0553*L;    %Cost based on Gunduz & Sahin No transmission and 100y flood Q
    %CostORNL = 1.2 * 110168 * DH^-0.35*(P*1e-3).^-0.3;   %Initial capital cost ORNL formula ($/kW)
    TurIrena{j} = 1e6 *(1.1943* (Pror*1e-6).^0.7634) * IR2005;  %Turbine investment costs ($) IRENA hydropower, empirical R2=0.94
    %     Pelton2jV = (1105.2*QQ.^-0.5106).*(P*1e-3)*ER;        %Pelton turbine Veileder p156 2-jetter high head Q<10 h=1000
    %     Pelton6jV = (1559.6*QQ.^-0.5179).*(P*1e-3)*ER;        %Pelton turbine Veileder p156 6-jetter high head 10>Q<35 h=1000 ($)
    %     FrancisV = (1439.1*QQ.^-0.3143).*(P*1e-3)*ER;         %Francis turbine Veileder p157 0>Q<160 h=300 ($)
    %     KaplanV = (11730.5*QQ.^-0.2953).*(P*1e-3)*ER;         %Kaplan turbine Veileder p158 0>Q<400 h=15 ($)
    %     SpecTurcost = TurIrena./(P*1e-3);                     %Specific investment costs ($)
    
    %Tunnel blasted costs based on Veileder report p67
    CorrMP = 0.0054*(L*1e-3)^2 - 0.0039*(L*1e-3) + 0.9671; %lenght multiplier correlation
    PCost = (219.99*Atunnel + 13658) * ER;   %Total price in ($/m2)
    PCostL = PCost * CorrMP;           %Cost times lenght multiplier ($/m)
    TPCostTunnel = PCostL*(L-Hgross);          %Total costs ($)
    %Steel penstock in tunnel costs based on Veileder report p94
    PenstockCost = (6*D + 9.4)*1e3*ER*Hgross; % $ Surface Penstock ($)
    TPCostL{j} = (TPCostTunnel + PenstockCost) * npipe(j);
    
    %Underground powerstation (depends on blasted volume) Veileder p101
    %Powerstation surface Veileder p106
    PSSnok = -0.0006*Qdesign.^2 + 0.67*Qdesign - 6.95;    % mil NOK/m (head 10-40m, 1 powerunit)
    PSSnok(PSSnok<20) = 20;                     % Not lower than 20mil NOK/m
    PSS(j) = PSSnok*ER*1e6;                        % $
    PSSkW = PSS(j)./(Pror*1e-3);                   % $/kW
    
    %Total electrical mechanical costs based on Veileder p145
    PEM{j} = (3.9142*(Pror*1e-6).^0.6622) * ER *1e6;       %Total electro-technical equipment $
    PEMkW = PEM{j}./(Pror*1e-3);                           %Total electro-technical equipment $/kW
    
    %Miscellaneous
    Misc{j} = (-38.795.*log(Qdesign) + 309.89).*(Pror*1e-3)*ER; %$
    
    %DisCost
    for i=1:numel(D)
        DisCostNOK(i) = Powerline_allocator(Dis,Pror(i)*1e-6); % NOK
        DisCost{j}(i) =DisCostNOK(i)*ER; % $2010
    end
    
    %% Collecting cost information
    AnnFac = interest/(1-((1+interest)^(-lifetime)));   % Annuity factor using interest (10%) and economic lifetime (40 years) = 10,23%
    
    AnnualCost = AnnFac * (TurIrena{j} + PEM{j} + PSS(j) + TPCostL{j} + Misc{j} + DisCost{j});
    COE = AnnualCost./PnetRoR;                           % Cost of electricity ($/kWh)
    
    AnTur = AnnFac * TurIrena{j};    %Turbine (annualized $)
    AnPipe = AnnFac * TPCostL{j};     %Pipe (annualized $)
    AnElec = AnnFac * PEM{j};         %Electro-technical equipment ($)
    AnPS = AnnFac * PSS(j);           %Powerstation($)
    AnMisc = AnnFac * Misc{j};        %Misc
    AnDis = AnnFac * DisCost{j};       %DisCost
    AnOM = AnnFac* ((TurIrena{j}+TPCostL{j}+PEM{j}+PSS(j)+Misc{j}+DisCost{j}) * OMshare);         %OM cost as fraction of total investment costs
    AnCtot = AnTur + AnPipe + AnElec + AnPS + AnMisc + AnOM + AnDis;
    AnOwner = AnCtot * Ownersrate;
    AnQuake = AnCtot * Quakerate;
    AnCtot = AnCtot + AnOwner + AnQuake;
    
    AnCtot_end{j} = AnCtot;
    
    COETur = AnTur./PnetRoR;
    COEOM = AnOM./PnetRoR;
    COEPipe = AnPipe./PnetRoR;
    COEElec = AnElec./PnetRoR;
    COEPS = AnPS./PnetRoR;
    COEMisc = AnMisc./PnetRoR;
    COEDis = AnDis./PnetRoR;
    COEOwner = AnOwner./PnetRoR;
    COEQuake = AnQuake./PnetRoR;
    COETot_np{j} = COETur+COEOM+COEPipe+COEElec+COEPS+COEMisc+COEDis+COEOwner+COEQuake+TransCost;
    
    %% Find optimal pipe diameter
    
    [COEminRoR_np(j) idx(j)]=min(COETot_np{j});
    
    Prormin(j) = Pror_end{j}(idx(j));
    hfmin(j) = hf_end{j}(idx(j));
    AnCtotmin(j) = AnCtot_end{j}(idx(j));
    
    PnetRoRmin_np(j) = PnetRoR(idx(j));
    %
    CostElements_np{j}(1)=COETur(idx(j));
    CostElements_np{j}(2)=COEOM(idx(j));
    CostElements_np{j}(3)=COEPipe(idx(j));
    CostElements_np{j}(4)=COEElec(idx(j));
    CostElements_np{j}(5)=COEPS(idx(j));
    CostElements_np{j}(6)=COEMisc(idx(j));
    CostElements_np{j}(7)=COEDis(idx(j));
    CostElements_np{j}(8)=COEOwner(idx(j));
    CostElements_np{j}(9)=COEQuake(idx(j));
    
    TurIrenaMin(j) = TurIrena{j}(idx(j));
    TPCostLMin(j) = TPCostL{j}(idx(j));
    PEMMin(j) = PEM{j}(idx(j));
    PSSMin(j) = PSS(j);
    MiscMin(j) = Misc{j}(idx(j));
    DisCostMin(j) = DisCost{j}(idx(j));
    
    if cost_constr==1
        
        if COEminRoR_np(j)>cost_lim; %Remove damsystem higher than x$/kWh
            COEminRoR_np(j)=NaN;
            PnetRoRmin_np(j)=NaN;
            CostElements_np{j}=NaN;
        else
        end
        
        %fprintf('no #%d, h #%d, in #%d: COE = %.2f $/kWh and Capacity = %.2f MW\n',k, heights, inletsp, COETotMx2, PMx2*1e-6)
        
    else
        %fprintf('no #%d, h #%d, in #%d: COE = %.2f $/kWh and Capacity = %.2f MW\n',k, heights, inletsp, COETotMx2, PMx2*1e-6)
        
    end
    
%     fprintf('Power is %.2f MWh.\n',PnetRoRmin_np(j)*1e-6);
%     fprintf('Lowest price is %.2f $/kWh.\n',COEminRoR_np(j));
%     fprintf('Optimal pipe diameter %.2f m.\n\n',D(idx(j)));
    
    %%
%     figure(1);
%     plot(D,COETot_np{j},'color',cc(j,:), 'DisplayName',['Npipe = ',num2str(npipe(j))],'LineWidth',1.2);
%     legend('-DynamicLegend')
%     title('Production price');
%     ylabel('$/kWh');
%     xlabel('Pipe diameter');
%     hold on
%     
%     figure(2);
%     bar(D,[COETur',COEOM', COEPipe',COEElec',COEPS', COEMisc', COEDis', COEOwner'],0.5,'stack')
%     legend('Turbine','OM','Pipe','Electro','PS','Misc','Dis','Owner')
%     legend('-DynamicLegend','location','Best')
%     title('COE stacked','FontWeight','bold');
%     ylabel('$/kWh');
%     xlabel('Pipe diameter (m)');
%     
%     %%
%     figure(3);
%     plot(D,PnetRoR*1e-6,'color',cc(j,:), 'DisplayName',['Npipe = ',num2str(npipe(j))],'LineWidth',1.2);
%     legend('-DynamicLegend')
%     title('Energy production');
%     ylabel('MWh');
%     xlabel('Pipe diameter');
%     hold on
    
end

[COEminRoR ind] = min(COEminRoR_np); % $/kWh
CostElements = CostElements_np{ind};
PnetRoRmin = PnetRoRmin_np(ind); %kWh
npipemin = npipe(ind); 
OptP = Prormin(ind); %W
OptInv = AnCtotmin(ind); % $ annualized
Opthf = hfmin(ind); % percentage of gross head lost due to friction loss (head loss)
%OptSpecCap = OptInv/(OptP*1e-3); % $/kW Annualized speccapcost
OptDP = idx(ind);

OptInv1 = TurIrenaMin(ind) + TPCostLMin(ind) + PEMMin(ind) + PSSMin(ind) ...
    + MiscMin(ind) + DisCostMin(ind); % $ absolute investments
OptInvOM = OMshare * OptInv1;
OptInvsOwners = (OptInv1 + OptInvOM) * Ownersrate;
OptInvsQuake = (OptInv1 + OptInvOM) * Quakerate;
OptInv2 = OptInv1 + OptInvOM + OptInvsOwners + OptInvsQuake;
 
OptSpecCap = OptInv2/(OptP*1e-3); % $/kW not annualized 

if COEminRoR<0.3; 
% fprintf('Number of pipes: %.0f.\n',npipe(ind));            
end

% fprintf('Power: %.2f MWh.\n',PnetRoRmin*1e-6);
% fprintf('Lowest price: %.2f $/kWh.\n',COEminRoR);
% fprintf('Number of pipes: %.0f.\n',npipe(ind));

%%

% (((OptInv2 * AnnFac)/(OptP*1e-3))  / (LF*8760*eta)) + TransCost

