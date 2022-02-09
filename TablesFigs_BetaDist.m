%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TablesFigs_BetaDist.m: Creates figures from Output tables
%% Dependencies: csv files in Data Files, SampleResPCs.mat, BilckeEtAlOutput.mat
%% ImportCovData, importfileBaseline, importfileCamp15, ImportResData73_BetaDist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dbstop if error
Folder = fileparts(pwd);
Folder = fullfile(Folder,'CleanPubCode/Figs');
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultAxesFontSize', 8)
set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultTextFontSize', 10)

%% Import Case and Resistance Data, Coverage data, and country R0 and reporting proportion
Baseline10Years73Total = importfileBaseline('DataFiles/Baseline10Years_73_TotalTyphoidNR.csv', [2, Inf]);
Baseline10Years73Fq = importfileBaseline('DataFiles/Baseline10Years_73_FQNSNR.csv', [2, Inf]);
Baseline10Years73MDR = importfileBaseline('DataFiles/Baseline10Years_73_MDRNR.csv', [2, Inf]);
Camp15Impact73Total = importfileCamp15('DataFiles/Camp15Impact_73_TotalTyphoidNR.csv', [2, Inf]);
Camp15Impact73Fq = importfileCamp15('DataFiles/Camp15Impact_73_FQNSNR.csv', [2, Inf]);
Camp15Impact73MDR = importfileCamp15('DataFiles/Camp15Impact_73_MDRNR.csv', [2, Inf]);

%Load in Beta distributions for countries
ImportResData73_BetaDist
%Load in Vaccine coverage
ImportCovData
%load in simulation results
load('SampleResPCs')
%Load in modeling estimates from previous study (Bilcke et al)
load('BilckeEtAlOutput.mat')
%Load in file with mean DR estimates for each country
load('DRtable.mat')

%% Calculate  proportion of Cases Averted for FQNS and MDR
FqPCMean = 100*Camp15Impact73Fq.CasesAvertedMedian./(Baseline10Years73Fq.CasesMedian);
FqPCMin = 100*Camp15Impact73Fq.CasesAvertedMin./(Baseline10Years73Fq.CasesMin);
FqPCMax = 100*Camp15Impact73Fq.CasesAvertedMax./(Baseline10Years73Fq.CasesMax);

FqOutput.PCMean=FqPCMean;
FqOutput.PCMin=FqPCMin;
FqOutput.PCMax=FqPCMax;
FqOutput.ISO3=CountryRegion.ISO3;

MDRPCMean = 100*Camp15Impact73MDR.CasesAvertedMedian./(Baseline10Years73MDR.CasesMedian);
MDRPCMin = 100*Camp15Impact73MDR.CasesAvertedMin./(Baseline10Years73MDR.CasesMin);
MDRPCMax = 100*Camp15Impact73MDR.CasesAvertedMax./(Baseline10Years73MDR.CasesMax);

MDROutput.PCMean=MDRPCMean;
MDROutput.PCMin=MDRPCMin;
MDROutput.PCMax=MDRPCMax;
MDROutput.ISO3=CountryRegion.ISO3;

FqPCMeanDeaths = 100*Camp15Impact73Fq.DeathsAvertedMedian./(Baseline10Years73Fq.DeathsMedian);
MDRPCMeanDeaths = 100*Camp15Impact73MDR.DeathsAvertedMedian./(Baseline10Years73MDR.DeathsMedian);
FqPCMeanDALYs = 100*Camp15Impact73Fq.DALYsAvertedMedian./(Baseline10Years73Fq.DALYsMedian);
MDRPCMeanDALYs = 100*Camp15Impact73MDR.DALYsAvertedMedian./(Baseline10Years73MDR.DALYsMedian);

%% % Marina's code from CEA adapted to include AMR data for maps
% set up folder name to save figures in correct folder

% Load in world data
load('world_map.mat') % the shp file, not a raster dataset.


%%
% Match initial DR and Impact with ISO country codes because country names are spelled
% differently
%

for i=1:size(world,1)
if sum(strcmp(cellstr(DRtable.ISO3), world(i,1).ISO3))==1
    world(i,1).Fq = .01*DRtable.FqMean(strcmp(cellstr(DRtable.ISO3), world(i,1).ISO3));
    world(i,1).MDR = .01*DRtable.MDRMean(strcmp(cellstr(DRtable.ISO3), world(i,1).ISO3));
    %if world(i,1).Fq ==0
     %   world(i,1).FqImp=0; %set impact to 0 rather than NaN if prevalence is 0
    %else
        world(i,1).FqImp = .01*FqOutput.PCMean(strcmp(cellstr(FqOutput.ISO3), world(i,1).ISO3));
    %end
    %if world(i,1).MDR==0
     %   world(i,1).MDRImp=0;
    %else
        world(i,1).MDRImp = .01*MDROutput.PCMean(strcmp(cellstr(MDROutput.ISO3), world(i,1).ISO3));
    %end
    %world(i,1).amr = countrycodes.amr(strcmp(countrycodes.ISO3, world(i,1).ISO3));
else 
    world(i,1).Fq = NaN; 
    world(i,1).MDR = NaN;
    world(i,1).FqImp = NaN; 
    world(i,1).MDRImp = NaN;    
end
end
%% Make new row for India bc plotting issue
India = world(98,:);
India.Lon=world(98,:).Lon(175:end);
India.Lat=world(98,:).Lat(175:end);

%% Make figure for Initial Resistant percentages
% Fq:
% Map

% MAP 1: Plot some spots on a map
figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,2,[1 3])
land = shaperead('landareas', 'UseGeoCoords', true);

% make a plain map of the world that serves as a background and plot some
% spots...
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,100])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
%set colormap
colormap(parula(73))
set(gca,'FontSize',9)

densityColorsFq = makesymbolspec('Polygon', {'Fq', [0 1], 'FaceColor', colormap}, {'Default', 'FaceColor', [0.8,0.8,0.8]});
axesm ('robinson', 'Frame', 'on', 'Grid', 'on')
geoshow(world, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsFq)
geoshow(India, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsFq)
c= colorbar;
ylabel(c,'Percent Resistant')

% Bar plots
subplot(4,2,2)
bar(DRtable.FqMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],DRtable.FqMean(1:37),DRtable.FqMean(1:37)-DRtable.FqMin(1:37),DRtable.FqMax(1:37)-DRtable.FqMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100],'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(4,2,4)
bar(DRtable.FqMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],DRtable.FqMean(38:73),DRtable.FqMean(38:73)-DRtable.FqMin(38:73),DRtable.FqMax(38:73)-DRtable.FqMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

% MDR
subplot(4,2,[5 7])
land = shaperead('landareas', 'UseGeoCoords', true);

% make a plain map of the world that serves as a background and plot some
% spots...
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,100])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
set(gca,'FontSize',9)

colormap(parula(73))
India = world(98,:);
India.Lon=world(98,:).Lon(175:end);
India.Lat=world(98,:).Lat(175:end);


densityColorsMDR = makesymbolspec('Polygon', {'MDR', [0 1], 'FaceColor', colormap}, {'Default', 'FaceColor', [0.8,0.8,0.8]});
axesm ('robinson', 'Frame', 'on', 'Grid', 'on')
geoshow(world, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsMDR)
geoshow(India, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsMDR)
c= colorbar;
ylabel(c,'Percent Resistant')

subplot(4,2,6)
bar(DRtable.MDRMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],DRtable.MDRMean(1:37),DRtable.MDRMean(1:37)-DRtable.MDRMin(1:37),DRtable.MDRMax(1:37)-DRtable.MDRMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(4,2,8)
bar(DRtable.MDRMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],DRtable.MDRMean(38:73),DRtable.MDRMean(38:73)-DRtable.MDRMin(38:73),DRtable.MDRMax(38:73)-DRtable.MDRMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100],'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)
%savefig(gcf,[Folder,'/Fig1_InitialAMR.fig'])
%saveas(gcf,[Folder,'/Fig1_InitialAMR'],'jpg')

%%
gtext('A)','FontSize',10,'FontWeight','Bold')
gtext('B)','FontSize',10,'FontWeight','Bold')
gtext('C)','FontSize',10,'FontWeight','Bold')
gtext('D)','FontSize',10,'FontWeight','Bold')

%% Make figure for Impact of vaccination on Fq Resistance
% Fq
% Map

% MAP 1: Plot some spots on a map
figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,2,[1 3])
land = shaperead('landareas', 'UseGeoCoords', true);

% make a plain map of the world that serves as a background and plot some
% spots...
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,100])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
%set colormap
colormap(autumn(1000))
set(gca,'FontSize',9)

densityColorsFqImp = makesymbolspec('Polygon', {'FqImp', [0 1], 'FaceColor', colormap}, {'Default', 'FaceColor', [0.8,0.8,0.8]});
axesm ('robinson', 'Frame', 'on', 'Grid', 'on')
geoshow(world, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsFqImp)
% India = world;
% world(98,:).Lon=world(98,:).Lon(175:end);
% world(98,:).Lat=world(98,:).Lat(175:end);
geoshow(India, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsFqImp)
c= colorbar;%('Limits',[0 1]); % for legend
%set(c,'Limits',[0 1])
ylabel(c,'Percent Averted')

% Bar chart
subplot(4,2,2)
bar(FqPCMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],FqPCMean(1:37),FqPCMean(1:37)-FqPCMin(1:37),FqPCMax(1:37)-FqPCMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100],'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(4,2,4)
bar(FqPCMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],FqPCMean(38:73),FqPCMean(38:73)-FqPCMin(38:73),FqPCMax(38:73)-FqPCMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100],'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

% MDR
% Map

% MAP 1: Plot some spots on a map
%figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,2,[5 7])
land = shaperead('landareas', 'UseGeoCoords', true);

% make a plain map of the world that serves as a background and plot some
% spots...
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,100])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
%set colormap
colormap(autumn(1000))
set(gca,'FontSize',9)

densityColorsMDRImp = makesymbolspec('Polygon', {'MDRImp', [0 1], 'FaceColor', colormap}, {'Default', 'FaceColor', [0.8,0.8,0.8]});
axesm ('robinson', 'Frame', 'on', 'Grid', 'on')
geoshow(world, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsMDRImp)
% India = world;
% world(98,:).Lon=world(98,:).Lon(175:end);
% world(98,:).Lat=world(98,:).Lat(175:end);
geoshow(India, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsMDRImp)
c= colorbar;%('Limits',[0 1]); % for legend
%set(c,'Limits',[0 1])
ylabel(c,'Percent Averted')

subplot(4,2,6)
bar(MDRPCMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],MDRPCMean(1:37),MDRPCMean(1:37)-MDRPCMin(1:37),MDRPCMax(1:37)-MDRPCMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100],'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(4,2,8)
bar(MDRPCMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],MDRPCMean(38:73),MDRPCMean(38:73)-MDRPCMin(38:73),MDRPCMax(38:73)-MDRPCMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100],'FontSize',9)
ytickformat(gca, 'percentage');
xtickangle(45)

%savefig(gcf,[Folder,'/Fig2_VaccImpactAMR.fig'])
%saveas(gcf,[Folder,'/Fig2_VaccImpactAMR'],'jpg')

%%
gtext('A)','FontSize',10,'FontWeight','Bold')
gtext('B)','FontSize',10,'FontWeight','Bold')
gtext('C)','FontSize',10,'FontWeight','Bold')
gtext('D)','FontSize',10,'FontWeight','Bold')


%% Make graphs of vaccine coverage vs cases averted, analyze with average age
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultTextFontSize', 12)

%Get mean coverage by country
CovTable = varfun(@nanmean,Gavicovunc,'InputVariables','Coverage','GroupingVariables','ISO3');

CountryRegion.Country = categorical(CountryRegion.Country);
CountryRegion.ISO3=categorical(CountryRegion.ISO3);
%add ISO to Camp15Impact
BL15Impact = outerjoin(Baseline10Years73Total,Camp15Impact73Total,'MergeKeys',true);
BL15Impact.Country=categorical(BL15Impact.country);
BL15Impact = outerjoin(BL15Impact,CountryRegion,'MergeKeys',true);

%merge coverage with impact
CovImpact = outerjoin(BL15Impact,CovTable,'MergeKeys',true);
CovImpact.PCAverted = 100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian;
for i = 1:size(CovImpact,1)
    CovImpact.AvgAgeInf(i) = avgage_country(find(cn ==find(categorical(iso.countryiso)==CovImpact.ISO3(i))));
    CovImpact.AgeDist(i) = agedist_country(find(cn ==find(categorical(iso.countryiso)==CovImpact.ISO3(i))));

end

CovImpact.AgeDist = categorical(CovImpact.AgeDist);
CovImpact.nanmean_Coverage=100*CovImpact.nanmean_Coverage
CovImpact_mdl = fitlm(CovImpact,'PCAverted~nanmean_Coverage+AvgAgeInf+AgeDist');

CovImpact_mdl_agedist = fitlm(CovImpact,'PCAverted~nanmean_Coverage+AgeDist');
CovImpact_mdl_agedist_int = fitlm(CovImpact,'PCAverted~nanmean_Coverage+AgeDist+nanmean_Coverage*AgeDist');

CovImpact_mdl_step= step(CovImpact_mdl);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
plot(CovImpact.nanmean_Coverage,100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian,'.','MarkerSize',15)
ytickformat(gca,'percentage')
xtickformat(gca,'percentage')
set(gca,'YLim',[0 100],'XLim',[30 100])
xlabel('Routine Vaccination Coverage')
ylabel('Percentage of Cases Averted')
%title("Percentage of Cases Averted by Vaccination Coverage")
subplot(2,3,2)
plot(CovImpact.AvgAgeInf,100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian,'.','MarkerSize',15)
ytickformat(gca,'percentage')
set(gca,'YLim',[0 100])
xlabel('Average Age of Infection')
%ylabel('Percentage of Cases Averted')
%title("Percentage of Cases Averted by Average Age of Infection")
subplot(2,3,3)
plot(CovImpact.AgeDist,100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian,'.','MarkerSize',15)
ytickformat(gca,'percentage')
set(gca,'YLim',[0 100],'XTickLabel',{'Very Young','Young','Average'})
xlabel('Age Distribution Category')
%ylabel('Percentage of Cases Averted')
%title("Percentage of Cases Averted by Average Age of Infection")
% saveas(gcf,[Folder,'/VaccCovAgePCAverted.fig'])
% saveas(gcf,[Folder,'/VaccCovAgePCAverted'],'jpg')

agestr={'Very Young','Young','Average'};
%figure('units','normalized','outerposition',[0 0 1 .5])
for i = 1:3
    subplot(2,3,i+3)
    plot(CovImpact.nanmean_Coverage(find(CovImpact.AgeDist==num2str(i))),100*CovImpact.CasesAvertedMedian(find(CovImpact.AgeDist==num2str(i)))./CovImpact.CasesMedian(find(CovImpact.AgeDist==num2str(i))),'.','MarkerSize',15)
    ytickformat(gca,'percentage')
    xtickformat(gca,'percentage')
    set(gca,'YLim',[0 100],'XLim',[30 100])
    xlabel('Routine Vaccination Coverage')
    if i==1
        ylabel('Percentage of Cases Averted')
    end
    title([agestr(i)])
end
sgtitle('Percentage of Cases Averted by Various Predictors','FontSize',24,'FontWeight','Bold')
gtext('A)','FontSize',16,'FontWeight','Bold')
gtext('B)','FontSize',16,'FontWeight','Bold')
savefig(gcf,[Folder,'/FigS5_CorrVaccImpact.fig'])
saveas(gcf,[Folder,'/FigS5_CorrVaccImpact'],'jpg')

%% Make plot of Difference in Resistance Percentage after intervention

%load in mat files for 10 years, 5 years, 20 years and make plot for each
% 10 years is already loaded: 
load('SampleResPCs')
figure('units','normalized','outerposition',[0 0 1 1])
subplot(3,1,1)
histogram(-100*DiffVec) %,'Normalization','probability')
set(gca,'XLim',[-20,80])
xtickformat(gca,'percentage')
ylabel('Number of Samples')
title('10 Years','FontSize',14)

load('SampleResPCsRes5Years')
subplot(3,1,2)
histogram(-100*DiffVec) %,'Normalization','probability')
set(gca,'XLim',[-20,80])
xtickformat(gca,'percentage')
ylabel('Number of Samples')
title('5 Years','FontSize',14)

load('SampleResPCsRes20Yrs')
subplot(3,1,3)
histogram(-100*DiffVec) %,'Normalization','probability')
set(gca,'XLim',[-20,80])
xtickformat(gca,'percentage')
ylabel('Number of Samples')
title('20 Years','FontSize',14)
xlabel('Percentage decrease in the relative prevalence of AMR')

gtext('A)','FontSize',16,'FontWeight','Bold')
gtext('B)','FontSize',16,'FontWeight','Bold')
gtext('C)','FontSize',16,'FontWeight','Bold')
sgtitle('Decrease in Proportion Resistant across 1200 Samples')
savefig(gcf,[Folder,'/FigS6_DiffResCamp1573.fig'])
saveas(gcf,[Folder,'/FigS6_DiffResCamp1573'],'jpg')


%% Mean percentage of cases averted (95% prediction interval)

%create matrix and write to excel file
%first row: FQ, Cases, Deaths Dalys
pcTable=[nanmean(FqOutput.PCMean) prctile(FqOutput.PCMean,2.5) prctile(FqOutput.PCMean,97.5)...
    nanmean(FqPCMeanDeaths) prctile(FqPCMeanDeaths,2.5) prctile(FqPCMeanDeaths,97.5)...
    nanmean(FqPCMeanDALYs) prctile(FqPCMeanDALYs,2.5) prctile(FqPCMeanDALYs,97.5);...
    nanmean(MDROutput.PCMean) prctile(MDROutput.PCMean,2.5) prctile(MDROutput.PCMean,97.5)...
    nanmean(MDRPCMeanDeaths) prctile(MDRPCMeanDeaths,2.5) prctile(MDRPCMeanDeaths,97.5)...
    nanmean(MDRPCMeanDALYs) prctile(MDRPCMeanDALYs,2.5) prctile(MDRPCMeanDALYs,97.5)];


