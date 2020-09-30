function [COETotminR OptPnet OptP OptDH OptDL CostElements OptInv Pmin OptPop OptLV OptSpecCap] = RDcostmodel(dfldamRD,dfhdamRD,dfPopCost,Qoutlets_design,Qoutlets_design_LF,Dis,LandValue,RDDepth,nbasin,outlet,cost_constr,cost_lim,Quakerate)

% b=k;
% dfldamRD=dfldamRD{b};
% dfhdamRD=single(dfhdamRD{b});
% dfPopCost=PopCost{b};
% Qoutlets_design=Qoutlets_design(b);
% Qoutlets_design_LF=Qoutlets_design_LF(b);
% Dis=DisOutlet(b);
% LandValue=LandValueRDlake{b};
% RDDepth=RDDepth(b);
% nbasin;
% outlet=b;
% cost_constr;
% cost_lim;
% Quakerate=0;

% Veileder 2012 Small hydropower http://publikasjoner.nve.no/veileder/2012/veileder2012_02.pdf

%% Variables from outer space
Qtot    = Qoutlets_design;      %Q_design based on monthly Q (4th highest month)
LF      = Qoutlets_design_LF;   %Loadfactor based on monthly pattern and Qdesign (30% exceedence)
DH      = dfhdamRD;             %Dam height (m)
DL      = dfldamRD;             %Dam lenght (m)
PopCost = dfPopCost;            %Lake cost (GDP * Pop)

%% constants
Hgross      = DH-RDDepth; %Gross head: Dam height minus riverdepth to determine actual pressurized head

for i=1:numel(DH)
    if Hgross(i)<3; Hgross(i)=0; end; % Minimum dam height of 3m
end

D           = 3;               %Pipe diameter (meters)
i           = 0;
cc          = hsv(10);         %Color generator
g           = 9.8;             %Gravitational acceleration (m/s2)
rho         = 1000;            %Density of water (kg/m3)
mu          = 0.001;           %Fluid viscosity of water (N-s/m2))
e           = 0.2;             %Roughness constant (m)
eta         = 0.7;             %Water to wire efficiency (turbine losses, pipe friction losses)
lifetime    = 40;              %Years
p           = 0.05;            %Price per kWh ($/kWh)
interest    = 0.1;             %Interest on capital
DW          = 15;              %Dam width (m)
OMshare     = 0.0225;          %As share of total investments
ER          = 0.1725;          %Exchange rate NOK to US$ 01/01/2010
IR2002      = 1.24;            %Inflation rate conversion to 2010 USD: World Bank Real effective exchange rate index (2010 = 100)
IR2005      = 1.09;            %Inflation rate conversion to 2010 USD: http://data.worldbank.org/indicator/PX.REX.REER?page=2
IR2013      = 0.99;            %Inflation rate conversion to 2010 USD: http://data.worldbank.org/indicator/PX.REX.REER
TransCost   = 0.0034 * IR2013; %$/kWh Fixed transmission cost per delivered kWh: http://www.eia.gov/forecasts/aeo/tables_ref.cfm Table 8 EIA ANNUAL ENERGY OUTLOOK 2015
Ownersrate  = 0.25;            %Owners cost due to lead times IRENA 2012 Hydropower

%% Physical calculations

QQ=Qtot;

PTgross = (rho*g*Qtot*1e-9*8760)*Hgross;      %Theoretical gross potential
PTnet = rho*g*Qtot*(Hgross)*8760*1e-3; % Theoretical net potential
P = eta*rho*g*Hgross.*QQ;                 % Capacity based on calculated Q (W)
Pnet = (P*1e-3*8760).*LF;              % Yearly production with LF (kWh)

%Different cost metrics based on capacity, Q, head and dam height
%CostGunduz = -0.0254 + 1.81*QQ + 0.0848*(Hgross-Ahf) + 0.0553*L;    %Cost based on Gunduz & Sahin No transmission and 100y flood Q
%CostORNL = 1.2 * 110168 * DH^-0.35*(P*1e-3).^-0.3;                  %Initial capital cost ORNL formula ($/kW)
TurIrena = 1e6 *(1.1943* (P*1e-6).^0.7634) * IR2005;                          %Turbine investment costs ($) IRENA hydropower, empirical R2=0.94
Pelton2jV = 1105.2*QQ.^-0.5106;                                       %Pelton turbine Veileder p156 2-jetter high head Q<10 h=1000
Pelton6jV = (1559.6*QQ.^-0.5179).*(P*1e-3) ;                          %Pelton turbine Veileder p156 6-jetter high head 10>Q<35 h=1000 ($/kW)
FrancisV = (1439.1*QQ.^-0.3143).*(P*1e-3);                                        %Francis turbine Veileder p157 0>Q<160 h=300 ($/kW)
KaplanV = 11730.5*QQ.^-0.2953;                                        %Kaplan turbine Veileder p158 0>Q<400 h=15 ($/kW)
SpecTurcost = TurIrena./(P*1e-3);                                    %Specific investment costs ($/kW)

%Dam cost (RCC) based on Veileder report p62
%DV = DL * 0.5*(DW*DH); %Dam volume (m3) 6e6 m3 (assumed it is a triangle)
DCn = 0.72*DH.^1.8;     %Price (1000 NOK/m) RCC >1M m3
%DCn = 1.69*DH.^1.68;   %Price (1000 NOK/m) RCC >0.1M m3
TDC = DCn.*DL*1e3*ER;   %Price ($)
TDC(TDC<20e6) =20e6;    %Minimum startup cost
TDamkW = TDC./(P*1e-3); %Dam cost per kW ($/kW)

%% Underground powerstation (depends on blasted volume) Veileder p101
%Powerstation surface Veileder p106
PSSnok = -0.0006*QQ.^2 + 0.67.*QQ - 6.95;    %mil NOK/m (head 10-40m, 1 powerunit)
PSSnok(PSSnok>180) = 180;                   % Not higher than 180mil NOK, levels off 
PSSnok(PSSnok<20) = 20;                     % Not lower than 20mil NOK
PSS = PSSnok*ER*1e6;                        % $
PSSkW = PSS./(P*1e-3);                      % $/kW

%% Total electrical mechanical costs based on Veileder p147 (2 power units)
PEM = (7.5748*(P*1e-6).^0.65618) * ER *1e6;       %Total electro-technical equipment $
PEMkW = PEM./(P*1e-3);                           %Total electro-technical equipment $/kW

%% Penstock cost Veileder p94 (surface penstock and steel pipes)
PENnok = (6*D+9.5)*DH.*(DL/65); %1000 NOK/m Veileder Eq. x Dam height x Number of pipes. Three Gorges has a pipe every 65m.
PEN = PENnok*ER*1e3;            %$
PENkW = PEN./(P*1e-3);          %$/kW

%% Fish pasage mitigation cost DOEwater p32: title: Estimation of Economic Parameters of U.S. Hydropower Resources
FishCost = 1.3e6 *(P*1e-6).^0.56 * IR2002;

%% Miscellaneous
Misc = (-38.795.*log(QQ) + 309.89).*(P*1e-3)*ER; % $

%% DisCost
DisCostNOK = Powerline_allocator(Dis,P*1e-6); % NOK
DisCost =DisCostNOK*ER; % $2010

%%
AnnFac = interest/(1-((1+interest)^(-lifetime)));   % Annuity factor using interest (10%) and economic lifetime (40 years) = 10,23%

AnnualCost = AnnFac * (TurIrena + PEM + PSS + TDC + PEN + FishCost + Misc + PopCost + DisCost + LandValue);
COE = AnnualCost./Pnet;                           % Cost of electricity ($/kWh)

%Collecting cost information
AnTurR = AnnFac * TurIrena;     %Turbine (annualized $)
AnDamR = AnnFac * TDC;          %Dam (annualized $)
AnPenR = AnnFac * PEN;          %Penstock (annualized $)
AnElecR = AnnFac * PEM;         %Electro-technical equipment ($)
AnPSR = AnnFac * PSS;           %Powerstation($)
AnFish = AnnFac * FishCost;     %Fish COst ($)
AnMisc = AnnFac * Misc;         %misc cost ($)
AnPopR = AnnFac * PopCost;      %Pop Cost ($)
AnDis = ((AnnFac * DisCost)*Misc)./Misc; %Distance2dem cost
AnLandVal = AnnFac * LandValue;  %Land value Cost ($)
AnOMR = AnnFac* ((TurIrena+TDC+PEN+PEM+PSS+FishCost+Misc+PopCost+DisCost+LandValue) * OMshare);         %OM cost as fraction of total investment costs
AnCtotR = AnTurR + AnDamR + AnPenR + AnElecR + AnPSR + AnFish + AnMisc + AnPopR + AnDis + AnLandVal + AnOMR;
AnCOwner = AnCtotR * Ownersrate;
AnCQuake = AnCtotR * Quakerate;
AnCtotR = AnCtotR + AnCOwner + AnCQuake;

COETurR = AnTurR./Pnet;
COEDamR = AnDamR./Pnet;
COEOMR = AnOMR./Pnet;
COEPenR = AnPenR./Pnet;
COEElecR = AnElecR./Pnet;
COEPSR = AnPSR./Pnet;
COEFish = AnFish./Pnet;
COEMisc = AnMisc./Pnet;
COEPop = AnPopR./Pnet;
COEDis = AnDis./Pnet;
COELandVal = AnLandVal./Pnet;
COEOwner = AnCOwner./Pnet;
COEQuake = AnCQuake./Pnet;
COETotR = COETurR+COEDamR+COEPenR+COEElecR+COEPSR+COEFish+COEMisc+COEPop+COEDis+COELandVal+COEOMR+COEOwner+COEQuake+TransCost;

%%
[c,b] =min(COETotR);
COETotminR = c;
OptDH = DH(b);
OptDL = DL(b);
OptP = P(b)*1e-6; %MW
OptPnet = Pnet(b); % kWh
OptInvAnn = AnCtotR(b); % $ Annualized investments
Pmin = P(b); % W
OptPop = PopCost(b); % $
OptLV = LandValue(b); % $

OptInv1 = TurIrena(b) + TDC(b) + PEN(b) + PEM(b) + PSS ...
    + FishCost(b) + Misc(b) + PopCost(b) + DisCost + LandValue(b); % $ absolute investments
OptInvOM = OMshare * OptInv1;
OptInvsOwners = (OptInv1 + OptInvOM) * Ownersrate;
OptInvsQuake = (OptInv1 + OptInvOM) * Quakerate;
OptInv = OptInv1 + OptInvOM + OptInvsOwners + OptInvsQuake;

OptSpecCap = OptInv/(Pmin*1e-3); % $/kW

CostElements{1}=COETurR;
CostElements{2}=COEDamR;
CostElements{3}=COEOMR;
CostElements{4}=COEPenR;
CostElements{5}=COEElecR;
CostElements{6}=COEPSR;
CostElements{7}=COEFish;
CostElements{8}=COEMisc;
CostElements{9}=COEPop;
CostElements{10}=COEDis;
CostElements{11}=COELandVal;
CostElements{12}=COEOwner;
CostElements{13}=COEQuake;
CostElements{14}=TransCost;

if cost_constr==1
    if COETotminR>cost_lim
        COETotminR=NaN; %Select river dam lower than x$/kWh
        OptPnet=NaN;
    else
    end
%     
%     h=figure(1);clf;
%     bar(1:numel(DH),[COETurR' COEDamR' COEOMR' COEPenR' COEElecR' COEPSR' COEFish' COEMisc' COEPop' COEDis' COELandVal' COEOwner'],0.5,'stack');
%     legend('Turbine', 'Dam', 'OM','Penstock','Electro','Powerstation','Fish','Misc','Pop','Distance','Land','Owners')
%     title('COE stacked');
%     ylabel('$/kWh');
%     xlabel('Dam height (m)');
%     axis([0,300,0,0.3])
%     
%     fprintf('Cost_constrain=on %d Outlet #%d, Cost = %.2f $/kWh, Cap = %.2f MW, DH = %.0f m\n',nbasin,outlet,COETotminR,OptP,OptDH)
%     
else
%     fprintf('Cost_constrain=off %d Outlet #%d, Cost = %.2f $/kWh, Cap = %.2f MW, DH = %.0f m\n',nbasin,outlet,COETotminR,OptP,OptDH)
    
end
% if COETotminR>2
%     COETotminR=NaN; %Select river dam lower than 0.5$/kWh
%     OptPnet=NaN;

% elseif COETotminR<0.2
%     h=figure('Visible','off');
%     bar(1:numel(DH),[COETurR' COEDamR' COEOMR' COEPenR' COEElecR' COEPSR' COEFish' COEMisc' COELandVal' COEDis'],0.5,'stack');
%     legend('Turbine', 'Dam', 'OM','Penstock','Electro','Powerstation','Fish','Misc','Lake','Distance')
%     title('COE stacked');
%     ylabel('$/kWh');
%     xlabel('Dam height (m)');
%     axis([0,inf,0,1])
%
%     pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
%     figfile = fullfile(pathname, sprintf('RDam_costperkWh_b%d_outlet%d.png',nbasin,outlet));
%     saveas(h,figfile);

%     hh=figure('Visible','off');
%     bar(1:numel(DH),[AnTurR' AnDamR' AnOMR' AnPenR' AnElecR' AnLandVal'],0.5,'stack');
%     legend('Turbine', 'Dam', 'OM','Penstock','Electro','Lake')
%     title('Investments stacked');
%     ylabel('Annualized investments');
%     xlabel('Dam width *90m');
%
%     pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
%     figfile = fullfile(pathname, sprintf('RDam_Invest_b%d_outlet%d.png',nbasin,outlet));
%     saveas(hh,figfile);

%     fprintf('Cost of outlet #%d = %.2f $/kWh\n',outlet, COETotminR)
%     fprintf('Capacity of outlet #%d = %.2f MW\n',outlet, OptP)
%     fprintf('DH of outlet #%d = %.0f m\n',outlet, OptDH)
%     fprintf('DL of outlet #%d = %.0f m\n',outlet, OptDL)
%     fprintf('Q of outlet #%d = %.2f m3/s \n',outlet, Qtot)
%     fprintf('LF of outlet #%d = %.2f \n\n',outlet, LF)
%     fprintf('Pnet of outlet #%d = %.0f MWh\n',outlet, OptPnet*1e-3)
% else
%
% end;

% COETotminRIdx = find(isnan(COETotminR)); %Find NaN index
% OptPnet(COETotminRIdx)=NaN;           %Remove Pnet with COE=NaN


%                 clf(figure(1),'reset');
%                 clf(figure(2),'reset');
%                 clf(figure(3),'reset');
%                 clf(figure(4),'reset');
%
%                 figure(1);
%                 plot(DH,COETot);
%                 title('COETot');
%                 xlabel('Dam width (m)');
%                 ylabel('$/kWh');
%                 %axis([500,1500,0,0.02]);
%
%                 figure(2)
%                 plot(DH,TDamkW);
%                 title('Dam cost');
%                 xlabel('Dam height (m)');
%                 ylabel('$/kW');
%
%
%                 figure(3)
%                 plot(DH,DL);
%                 title('Dam height and dam width');
%                 xlabel('Dam height (m)');
%                 ylabel('Dam width (m)');
%                 %%
% h=figure(3);
% bar(1:20,[COETurR' COEDamR' COEOMR' COEPenR' COEElecR' COELake'],0.5,'stack');
% legend('Turbine', 'Dam', 'OM','Penstock','Electro','Lake')
% title('COE stacked');
% ylabel('$/kWh');
% xlabel('Dam width (*90m)');
%
% pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
% figfile = fullfile(pathname, sprintf('RDam_costperkWh_b%d_outlet%d.png',nbasin,outlet));
% saveas(h,figfile);
%
% %%
% hh=figure(4);
% bar(1:20,[AnTurR' AnDamR' AnOMR' AnPenR' AnElecR' AnLakeR'],0.5,'stack');
% legend('Turbine', 'Dam', 'OM','Penstock','Electro','Lake')
% title('Investments stacked');
% ylabel('Annualized investments');
% xlabel('Dam width *90m');
%
% pathname = fileparts('Y:\ontwapps\Timer\Users\David\Pojects\Hydropower\Model\NAM\output\RD\');
% figfile = fullfile(pathname, sprintf('RDam_Invest_b%d_outlet%d.png',nbasin,outlet));
% saveas(hh,figfile);

% end
%% Formula to re-calculate starting with SpecCapCost

% (((OptInv * AnnFac)/(P(b)*1e-3))  / (LF*8760*eta)) + TransCost
