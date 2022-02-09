clear
dbstop if error
Folder = fileparts(pwd);
Folder = fullfile(Folder,'Slides','Figs');
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultAxesFontSize', 8)
set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultTextFontSize', 10)

%% Import Case and Resistance Data, Coverage data, and country R0 and reporting proportion
Baseline10Years73Total = importfileBaseline("Baseline10_73_Total.csv", [2, Inf]);
%Baseline10Years73Sens = importfileBaseline("Baseline10_73_Sensitive.csv", [2, Inf]);
Baseline10Years73Fq = importfileBaseline("Baseline10_73_FQNS.csv", [2, Inf]);
Baseline10Years73MDR = importfileBaseline("Baseline10_73_MDR.csv", [2, Inf]);
Camp15Impact73Total = importfileCamp15("Camp15Impact_73_Total.csv", [2, Inf]);
%Camp15Impact73Sens = importfileCamp15("Camp15Impact_73_Sensitive.csv", [2, Inf]);
Camp15Impact73Fq = importfileCamp15("Camp15Impact_73_FQNS.csv", [2, Inf]);
Camp15Impact73MDR = importfileCamp15("Camp15Impact_73_MDR.csv", [2, Inf]);


ImportResData73_BetaDist
ImportCovData
%load('SampleResPCs20200403T155735')
%load('SampleResPCs20200514T104422')
load('SampleResPCs20200606T212037')
load('MarinaFiles/MarinaOutput.mat')
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

%% % Marina's code from CEA adapted to include AMR data for maps
% set up folder name to save figures in correct folder

% Load in world data
load('Maps/world_map.mat') % the shp file, not a raster dataset.


%%
% Match initial DR and Impact with ISO country codes because country names are spelled
% differently
%

for i=1:size(world,1)
if sum(strcmp(cellstr(DRtable.ISO3), world(i,1).ISO3))==1
    world(i,1).Fq = .01*DRtable.FqMean(strcmp(cellstr(DRtable.ISO3), world(i,1).ISO3));
    world(i,1).MDR = .01*DRtable.MDRMean(strcmp(cellstr(DRtable.ISO3), world(i,1).ISO3));
    if world(i,1).Fq ==0
        world(i,1).FqImp=0; %set impact to rather than NaN if prevalence is 0
    else
        world(i,1).FqImp = .01*FqOutput.PCMean(strcmp(cellstr(FqOutput.ISO3), world(i,1).ISO3));
    end
    if world(i,1).MDR==0
        world(i,1).MDRImp=0;
    else
        world(i,1).MDRImp = .01*MDROutput.PCMean(strcmp(cellstr(MDROutput.ISO3), world(i,1).ISO3));
    end
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

densityColorsFq = makesymbolspec('Polygon', {'Fq', [0 1], 'FaceColor', colormap}, {'Default', 'FaceColor', [0.8,0.8,0.8]});
axesm ('robinson', 'Frame', 'on', 'Grid', 'on')
geoshow(world, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsFq)
% India = world;
% world(98,:).Lon=world(98,:).Lon(175:end);
% world(98,:).Lat=world(98,:).Lat(175:end);
geoshow(India, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsFq)
c= colorbar;%('Limits',[0 1]); % for legend
%set(c,'Limits',[0 1])
ylabel(c,'Percent Resistant')
gtext('A)','FontSize',16,'FontWeight','Bold')

% Bar plots
subplot(4,2,2)
bar(DRtable.FqMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],DRtable.FqMean(1:37),DRtable.FqMean(1:37)-DRtable.FqMin(1:37),DRtable.FqMax(1:37)-DRtable.FqMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100])
ytickformat(gca, 'percentage');
gtext('B)','FontSize',16,'FontWeight','Bold')
xtickangle(45)

subplot(4,2,4)
bar(DRtable.FqMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],DRtable.FqMean(38:73),DRtable.FqMean(38:73)-DRtable.FqMin(38:73),DRtable.FqMax(38:73)-DRtable.FqMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73))
ytickformat(gca, 'percentage');
xtickangle(45)
%sgtitle('B)','FontSize',24,'FontWeight','Bold')

%saveas(gcf,[Folder,'/MapFqPc73'],'jpg')

% MDR
%figure('units','normalized','outerposition',[0 0 1 1])
subplot(4,2,[5 7])
land = shaperead('landareas', 'UseGeoCoords', true);

% make a plain map of the world that serves as a background and plot some
% spots...
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,100])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])

colormap(parula(73))
India = world(98,:);
India.Lon=world(98,:).Lon(175:end);
India.Lat=world(98,:).Lat(175:end);

% 0-1 might be a crazy range if the data is mostly 0.1-0.2, so change that
% here

densityColorsMDR = makesymbolspec('Polygon', {'MDR', [0 1], 'FaceColor', colormap}, {'Default', 'FaceColor', [0.8,0.8,0.8]});
axesm ('robinson', 'Frame', 'on', 'Grid', 'on')
geoshow(world, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsMDR)
% India = world;
% world(98,:).Lon=world(98,:).Lon(175:end);
% world(98,:).Lat=world(98,:).Lat(175:end);
geoshow(India, 'DisplayType', 'polygon', 'SymbolSpec', densityColorsMDR)
c= colorbar;%('Limits',[0 1]); % for legend
%set(c,'Limits',[0 1])
ylabel(c,'Percent Resistant')
gtext('C)','FontSize',16,'FontWeight','Bold')

subplot(4,2,6)
bar(DRtable.MDRMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],DRtable.MDRMean(1:37),DRtable.MDRMean(1:37)-DRtable.MDRMin(1:37),DRtable.MDRMax(1:37)-DRtable.MDRMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37))
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(4,2,8)
bar(DRtable.MDRMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],DRtable.MDRMean(38:73),DRtable.MDRMean(38:73)-DRtable.MDRMin(38:73),DRtable.MDRMax(38:73)-DRtable.MDRMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100])
ytickformat(gca, 'percentage');
xtickangle(45)
gtext('D)','FontSize',16,'FontWeight','Bold')
%title("Initial Percent of Cases Resistant: MDR")
savefig(gcf,[Folder,'/Fig1_InitialAMR.fig'])
saveas(gcf,[Folder,'/Fig1_InitialAMR'],'jpg')

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
gtext('A)','FontSize',16,'FontWeight','Bold')

% Bar chart
subplot(4,2,2)
bar(FqPCMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],FqPCMean(1:37),FqPCMean(1:37)-FqPCMin(1:37),FqPCMax(1:37)-FqPCMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100])
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(4,2,4)
bar(FqPCMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],FqPCMean(38:73),FqPCMean(38:73)-FqPCMin(38:73),FqPCMax(38:73)-FqPCMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100])
ytickformat(gca, 'percentage');
xtickangle(45)
gtext('B)','FontSize',16,'FontWeight','Bold')


% Make figure for Impact of vaccination on MDR Resistance

% MAP 1: Plot some spots on a map
figure('units','normalized','outerposition',[0 0 1 1])
land = shaperead('landareas', 'UseGeoCoords', true);

% make a plain map of the world that serves as a background and plot some
% spots...
ax = axesm ('robinson', 'Frame', 'on', 'Grid', 'off');
caxis([0,100])
geoshow(ax, land, 'FaceColor', [0.8 0.8 0.8])
%set colormap
colormap(autumn(1000))

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
title("Percent Multidrug Resistant Cases Averted")
saveas(gcf,[Folder,'/MapMDRImp73'],'jpg')



%% Make graphs of resistance prevalence by country, bar chart and map

%Fq Plot
%figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1)
% bar(DRtable.FqMean(1:37),'FaceColor',[.6 .6 .6])
% hold on
% er = errorbar([1:37],DRtable.FqMean(1:37),DRtable.FqMean(1:37)-DRtable.FqMin(1:37),DRtable.FqMax(1:37)-DRtable.FqMean(1:37));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100])
% ytickformat(gca, 'percentage');
% xtickangle(45)
% 
% subplot(2,1,2)
% bar(DRtable.FqMean(38:73),'FaceColor',[.6 .6 .6])
% hold on
% er = errorbar([1:36],DRtable.FqMean(38:73),DRtable.FqMean(38:73)-DRtable.FqMin(38:73),DRtable.FqMax(38:73)-DRtable.FqMean(38:73));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73))
% ytickformat(gca, 'percentage');
% xtickangle(45)
% sgtitle('B)','FontSize',24,'FontWeight','Bold')
%sgtitle('Percent of Cases Fluoroquinolone Resistant by Country','FontSize',24,'FontWeight','Bold')
%saveas(gcf,[Folder,'/CountryFqPc73'],'jpg')
%savefig(gcf,[Folder,'/CountryFqPc73.fig'])

% MDR Plot
%figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1)
% bar(DRtable.MDRMean(1:37),'FaceColor',[.6 .6 .6])
% hold on
% er = errorbar([1:37],DRtable.MDRMean(1:37),DRtable.MDRMean(1:37)-DRtable.MDRMin(1:37),DRtable.MDRMax(1:37)-DRtable.MDRMean(1:37));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37))
% ytickformat(gca, 'percentage');
% xtickangle(45)
% 
% subplot(2,1,2)
% bar(DRtable.MDRMean(38:73),'FaceColor',[.6 .6 .6])
% hold on
% er = errorbar([1:36],DRtable.MDRMean(38:73),DRtable.MDRMean(38:73)-DRtable.MDRMin(38:73),DRtable.MDRMax(38:73)-DRtable.MDRMean(38:73));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100])
% ytickformat(gca, 'percentage');
% xtickangle(45)
%sgtitle('D)','FontSize',24,'FontWeight','Bold')
%sgtitle('Percent of Cases Multi-drug Resistant by Country','FontSize',24,'FontWeight','Bold')
%savefig(gcf,[Folder,'/CountryMDRPc73.fig'])
%saveas(gcf,[Folder,'/CountryMDRPc73'],'jpg')

%% Plot proportion of Cases Averted for FQNS and MDR

figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1)
% bar(FqPCMean(1:37),'FaceColor',[.6 .6 .6])
% hold on
% er = errorbar([1:37],FqPCMean(1:37),FqPCMean(1:37)-FqPCMin(1:37),FqPCMax(1:37)-FqPCMean(1:37));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100])
% ytickformat(gca, 'percentage');
% xtickangle(45)
% 
% subplot(2,1,2)
% bar(FqPCMean(38:73),'FaceColor',[.6 .6 .6])
% hold on
% er = errorbar([1:36],FqPCMean(38:73),FqPCMean(38:73)-FqPCMin(38:73),FqPCMax(38:73)-FqPCMean(38:73));    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100])
% ytickformat(gca, 'percentage');
% xtickangle(45)
% sgtitle('B)','FontSize',24,'FontWeight','Bold')
%sgtitle('Percent of Fluoroquinolone-Resistant Cases Averted by Country','FontSize',24,'FontWeight','Bold')
savefig(gcf,[Folder,'/CountryFqPcAverted73.fig'])
%saveas(gcf,[Folder,'/CountryFqPcAverted73'],'jpg')

% MDR Plot
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
bar(MDRPCMean(1:37),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:37],MDRPCMean(1:37),MDRPCMean(1:37)-MDRPCMin(1:37),MDRPCMax(1:37)-MDRPCMean(1:37));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:37],'xticklabels',DRtable.Country(1:37),'YLim',[0 100])
ytickformat(gca, 'percentage');
xtickangle(45)

subplot(2,1,2)
bar(MDRPCMean(38:73),'FaceColor',[.6 .6 .6])
hold on
er = errorbar([1:36],MDRPCMean(38:73),MDRPCMean(38:73)-MDRPCMin(38:73),MDRPCMax(38:73)-MDRPCMean(38:73));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
set(gca,'XTick',[1:36],'xticklabels',DRtable.Country(38:73),'YLim',[0 100])
ytickformat(gca, 'percentage');
xtickangle(45)
sgtitle('D)','FontSize',24,'FontWeight','Bold')
%sgtitle('Percent of MDR Cases Averted by Country','FontSize',24,'FontWeight','Bold')
savefig(gcf,[Folder,'/CountryMDRPcAverted73.fig'])
%saveas(gcf,[Folder,'/CountryMDRPcAverted73'],'jpg')

%% Make graphs of vaccine coverage vs cases averted, analyze with average age

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

%CovImpact.AgeDist = ordinal(CovImpact.AgeDist);

CovImpact_mdl = fitlm(CovImpact,'PCAverted~nanmean_Coverage+AvgAgeInf+AgeDist');

CovImpact_mdl= step(CovImpact_mdl);
figure('units','normalized','outerposition',[0 0 1 .5])
subplot(1,3,1)
plot(100*CovImpact.nanmean_Coverage,100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian,'.','MarkerSize',15)
ytickformat(gca,'percentage')
xtickformat(gca,'percentage')
set(gca,'YLim',[0 100],'XLim',[30 100])
xlabel('Routine Vaccination Coverage')
ylabel('Percentage of Cases Averted')
%title("Percentage of Cases Averted by Vaccination Coverage")
subplot(1,3,2)
plot(CovImpact.AvgAgeInf,100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian,'.','MarkerSize',15)
ytickformat(gca,'percentage')
set(gca,'YLim',[0 100])
xlabel('Average Age of Infection')
ylabel('Percentage of Cases Averted')
%title("Percentage of Cases Averted by Average Age of Infection")
subplot(1,3,3)
plot(CovImpact.AgeDist,100*CovImpact.CasesAvertedMedian./CovImpact.CasesMedian,'.','MarkerSize',15)
ytickformat(gca,'percentage')
set(gca,'YLim',[0 100],'XTick',[1:3],'XTickLabel',{'Very Young','Young','Average'})
xlabel('Age Distribution Category')
ylabel('Percentage of Cases Averted')
%title("Percentage of Cases Averted by Average Age of Infection")
sgtitle('Percentage of Cases Averted by Various Predictors','FontSize',24,'FontWeight','Bold')
saveas(gcf,[Folder,'/VaccCovAgePCAverted'],'jpg')

agestr={'Very Young','Young','Average'};
figure('units','normalized','outerposition',[0 0 1 .5])
for i = 1:3
    subplot(1,3,i)
    plot(100*CovImpact.nanmean_Coverage(find(CovImpact.AgeDist==i)),100*CovImpact.CasesAvertedMedian(find(CovImpact.AgeDist==i))./CovImpact.CasesMedian(find(CovImpact.AgeDist==i)),'.','MarkerSize',15)
    ytickformat(gca,'percentage')
    xtickformat(gca,'percentage')
    set(gca,'YLim',[0 100],'XLim',[30 100])
    xlabel('Routine Vaccination Coverage')
    ylabel('Percentage of Cases Averted')
    title(['Age Distribution: ', agestr(i)])
end
sgtitle('Cases Averted vs Vaccination Coverage by Population Age Distribution','FontSize',24,'FontWeight','Bold')
savefig(gcf,[Folder,'/CorrVaccImpact.fig'])
%saveas(gcf,[Folder,'/CorrVaccImpact'],'jpg')

% [corrCov,pCov] =corr(CovImpact.nanmean_Coverage,CovImpact.MedianCasesAverted./CovImpact.MedianCases);
% text(80,70,['corr = ',num2str(corr(CovImpact.nanmean_Coverage,CovImpact.MedianCasesAverted./CovImpact.MedianCases))]);
%% Make plot of Difference in Resistance Percentage after intervention
% figure('units','normalized','outerposition',[0 0 1 1])
% histogram(-100*DiffVec) %,'Normalization','probability')
% xtickformat(gca,'percentage')
% title('Decrease in Proportion Resistant across 1200 Samples','FontSize',24)
% saveas(gcf,[Folder,'/DiffResCamp1573_5Years'],'jpg')
%% Display numbers for paper


%% Mean percentage of cases averted (95% prediction interval)

%create matrix and write to excel file
%first row: FQ, Cases, Deaths Dalys
pcTable=[nanmean(FqOutput.PCMean) prctile(FqOutput.PCMean,2.5) prctile(FqOutput.PCMean,97.5)...
    nanmean(FqPCMeanDeaths) prctile(FqPCMeanDeaths,2.5) prctile(FqPCMeanDeaths,97.5)...
    nanmean(FqPCMeanDALYs) prctile(FqPCMeanDALYs,2.5) prctile(FqPCMeanDALYs,97.5);...
    nanmean(MDROutput.PCMean) prctile(MDROutput.PCMean,2.5) prctile(MDROutput.PCMean,97.5)...
    nanmean(MDRPCMeanDeaths) prctile(MDRPCMeanDeaths,2.5) prctile(MDRPCMeanDeaths,97.5)...
    nanmean(MDRPCMeanDALYs) prctile(MDRPCMeanDALYs,2.5) prctile(MDRPCMeanDALYs,97.5)];


%% Write output to csv files and mat files in maps folder

% writetable(FqOutput,'FqOutput73.csv')
% writetable(MDROutput,'MDROutput73.csv')

%save Outputs table into maps folder
save('Maps/ImpactOutputs.mat','FqOutput','MDROutput')
