%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OutputTablesDR_ResDiff: outputs a mat file with mean, min and max of FQNS
% and MDR estimates for each country for Figures, outputs table with matrix
% of changes in resistance proportion for each simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dbstop if error
Folder = fileparts(pwd);
Folder = fullfile(Folder,'CleanPubCode/Figs');

%% Import Case and Resistance Data, Coverage data, and country R0 and reporting proportion

%Load in Resistance estimates for countries
ImportResData73_Means
%Load in Vaccine coverage
ImportCovData
%load in simulation results
load('SampleResPCs')
%Load in modeling estimates from previous study (Bilcke et al)
load('BilckeEtAlOutput.mat')

%% Collate Resistance data
%replace Kenya Fq with 0 temporarily and save MDR number
CountryDRData.Fq(find(CountryDRData.Country=="Kenya"))=0;
KenyaMDR = CountryDRData.MDR(find(CountryDRData.Country=="Kenya"));
%replace Guinea-Bissau and Senegal MDR with 0 temporarily and save Fq
%number
CountryDRData.MDR(find(CountryDRData.Country=="Guinea Bissau" | CountryDRData.Country=="Senegal"))=0;
GBFq=CountryDRData.Fq(find(CountryDRData.Country=="Guinea Bissau"));
SenegalFq=CountryDRData.Fq(find(CountryDRData.Country=="Senegal"));

%% Calculate mean, min, and max and put in table
DRtable1 = grpstats(CountryDRData,'Country',{'nanmean','nanmin','nanmax'});
DRtable1 = removevars(DRtable1,'GroupCount');
DRtable1.Properties.VariableNames= {'Country','FqMean','FqMin','FqMax','MDRMean','MDRMin','MDRMax'};

%% Make Table of countrynames
CountryTab=table(Country);
CountryTab.Country=categorical(Country);
CountryTab.ISO3 = categorical(ISO3);

%Outer join with DR table
DRtable=outerjoin(CountryTab,DRtable1,'MergeKeys',true);

%% Get regional DR estimates
% Change to categorical variables
CountryRegion.Country=categorical(CountryRegion.Country);
RegionDRData.Region=categorical(RegionDRData.Region);
%Get regional estimates
CountryRegionDR = join(CountryRegion,RegionDRData);
CountryRegionDR=removevars(CountryRegionDR,'Region');

%% Join regional estimates to main table
[idxA, idxB] = ismember(DRtable.Country, CountryRegionDR.Country);
DRtable(idxA,{'FqMean','FqMin','FqMax','MDRMean','MDRMin','MDRMax'}) = CountryRegionDR(idxB(idxA),{'FqMean','FqMin','FqMax','MDRMean','MDRMin','MDRMax'});
%Fill in missing values for Kenya, Guinea-Bissau and Senegal
DRtable(find(DRtable.Country=="Kenya"),{'MDRMean','MDRMin','MDRMax'})={KenyaMDR,.95*KenyaMDR,1.05*KenyaMDR};
DRtable(find(DRtable.Country=="Guinea Bissau"),{'FqMean','FqMin','FqMax'})={GBFq,GBFq,GBFq};
DRtable(find(DRtable.Country=="Senegal"),{'FqMean','FqMin','FqMax'})={SenegalFq,.95*SenegalFq,1.05*SenegalFq};

clear DRtable1 idxA idxB KenyaMDR SenegalFq GBFq

%save DR table 
save('DRtable.mat','DRtable')
%% Get Change in Proportion of Resistant Cases

% interpolate change -- choose midpoint of change in resistant infections
% for range to multiply

% get uniform distribution from min and max of resistance percentages
% multiply times median, min, max of cases 

% multiply times uniform distribution of change in Resistance min max midpoints
% get 2000 samples per country of difference in resistance calculated using
% R0 and repcountry with resDiff_mdl

resDiffsMat = zeros(length(Country),2000);
for i = 1:length(Country)
    
    % get 2000 samples per country of difference in resistance calculated using
    % R0 and repcountry with resDiff_mdl
    % get correct R0 and prop symptomaticfrom country: 
    R0tmp = R0country(find(cn ==find(categorical(iso.countryiso)==DRtable.ISO3(i))));
    symptmp = repcountry(find(cn ==find(categorical(iso.countryiso)==DRtable.ISO3(i))));
    % use model to predict 2000 samples of resistance difference
    resDiffs = .01*random(resDiff_mdl,repmat([R0tmp,symptmp],2000,1));
    resDiffsMat(i,:) = resDiffs;
        
end

resDiffsTable=array2table(resDiffsMat);
resDiffsTable= [CountryTab resDiffsTable];

%make histograms of high and low R0 change in resistance
% low, Somalia

figure('units','normalized','outerposition',[0 0 1 .5])
subplot(1,2,1)
histogram(100*table2array(resDiffsTable(60,3:end)),'Normalization','probability')
xtickformat(gca,'percentage')
title('Range of Differences in Proportion Resistant, Somalia, R_0 = 1.49') 
subplot(1,2,2)
histogram(100*table2array(resDiffsTable(7,3:end)),'Normalization','probability')
xtickformat(gca,'percentage')
title('Range of Differences in Proportion Resistant, Bhutan, R_0 = 10.4') 
saveas(gcf,[Folder,'/FigS8_R0ResDiffFigs'],'jpg')

writetable(resDiffsTable,'resDiffstable.csv')