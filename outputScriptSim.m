%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputScriptSim.m: Script for running simulations for different R0s,
% proportion symptomatic, vax coverage, and proportion resistant. Runs
% model, generates plots for Figs X Y Z.
% Dependencies: sir_amr_age.m, q_only.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
dbstop if error

%Set up figure fonts
set(0,'DefaultAxesFontName','Times New Roman')
set(0,'DefaultAxesFontSize', 18)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 18)


%Set up folder to save figures
Folder = fileparts(pwd);
Folder = fullfile(Folder,'CleanPubCode/Figs');

% Load file with parameter values
load('q_only')

%Set up R0 values to sample:
R0samp = [1.5 2.5 3.5 7 10.5];

%Set up proportion symptomatic vector (using 5th pctile, median, mean, 95%
%pctile from previous fits, roughly)
sympSamp = [.01 .05 .10 .5];

%Set up vaccination coverage matrix
CovMat = [ .7:.01:.79; .62:.02:.8; .86:.01:.95];

% Set up resistance parameters 
% r: rate of resistance acquisitions
% d1r d2r durations of infectiousness (symptomatic, asymptomatic)
% rr: relative infectiousness 
resParams = struct;
% Latin Hypercube sample of resistance parameters
nsamp=20; %number of samples
lhs_sample=lhsdesign(nsamp,4);
param_max =[ 5 .7 .7 3]; % 
param_min = [ .01 .075 .075 .3];
param_sample=param_min+lhs_sample.*(param_max-param_min);

for i = 1:nsamp
    resParams(i).r = param_sample(i,1);
    resParams(i).d1r = param_sample(i,2);
    resParams(i).d2r = param_sample(i,3);
    resParams(i).rr = param_sample(i,4);
end

% age distributions from original model-- use lmid for these simulations
low_inc = [0.160015881	0.141894822	0.124871077	0.108908019	0.092259324	0.076602459];
lmid_inc = [0.1061743	0.103271048	0.099353973	0.095006697	0.090143179	0.08454642]; 
umid_inc = [0.071794938	0.069233491	0.067531964	0.069770015	0.075625419	0.089103213]; 

dist_mat = [low_inc;lmid_inc;umid_inc];

% Set up birthrates for low, mid, upper mid income countries
birthrates = 1./[36.6, 23.6, 15.0];


%Pre-allocate output matrices
CasesAvertedMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),3);
ResCasesAvertedMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),3);
ResCarAvertedMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),3);
DeathsAvertedMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),3);
AbxAvertedMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),3);
PercentResCasesMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),4);
PercentResCasesMatDN = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),4);
PercentResCasesMatTrDN = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),4);
PercentSensCasesMat = zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),4);
ResCasesTime= zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),780,4);
SensCasesTime= zeros(size(R0samp,2),size(sympSamp,2),size(CovMat,1),size(resParams,2),780,4);
%% Run loops for different R0s, proportions symptomatic, vax coverages, initial resistance levels
for j = 1:size(R0samp,2)
    for k = 1:size(sympSamp,2)
        for l = 1:size(CovMat,1)
            for m = 1:size(resParams,2)
            

%simulation lengths 
t0=5200; %length of "burn-in period" (in weeks)
tvacc=260; %time of vaccine introduction


% Set distributions to Low Middle Income for simulations
dist_num = 2; 
dist = dist_mat(dist_num,:);

% strata age groups with young group breakdown for vax: 0-9 months,
% 9months-2 years, 2-5, 5-9, 10-14, 15-29, 30+
agedist = [3*dist(1)/20; dist(1)/4; dist(1)*3/5; dist(2); dist(3); sum(dist(4:6)); 1-sum(dist)];

B=birthrates(dist_num); 

%Set up parameters
q.mub=log(1+B)/52.18; % birthrates
q.mu = q.mub; 
q.N = 1e5; %pop size of 100,000
q.theta = sum(agedist(1:6)*q.t1s_y + agedist(7)*q.t1s_o); %Carriers
q.betak = R0samp(j).*(q.mub + q.d1s)./((1+q.d1s*q.theta*q.rc./q.mub)); %force of infection
q.rr_age1 = .37; %relative risk of youngest age group
q.rr_age2 = .68; %relative risk of second youngest age group
q.s = sympSamp(k); %proportion symptomatic
%resistance parameters
q.r = resParams(m).r;
q.d1r = resParams(m).d1r;
q.d2r = resParams(m).d2r;
q.rr =resParams(m).rr;



%vaccination
rout_cov1 = CovMat(l,:); %Coverage over 10 years
rout_cov = repmat(rout_cov1,52,1);
rout_cov = reshape(rout_cov,length(rout_cov1)*52,1);
tf=t0+tvacc+length(rout_cov); %total length of simulation
v=zeros(tf,1);

%campaign vaccination scenarios
camp_cov = .9;
camp_age = 1-(1-camp_cov)^.25;
mass_vac = zeros(tf,q.n_age);
%age groups
vacc_age = camp_age*[0 0 0 0 0 0 0; 0 1 1 1  0 0 0; 0 1 1 1 1 0 0];
camp_label={'No Vaccination','Routine Only','Campaign to 5 yrs','Campaign to 15 yrs'};

population =q.N*agedist;

%Loop through vaccination coverage scenarios:
% 1: No Vax, 2: Routine only, 
%3: routine + catch-up to 5, 4: routine +catch-up to 15
for c=1:4 
if c >1
    v=[zeros(t0+tvacc,1); rout_cov];
    mass_vac = [zeros(t0+tvacc,q.n_age);repmat(vacc_age(c-1,:),4,1);zeros(tf-(t0+tvacc+4),q.n_age)];
end
%State of the population
St0=zeros(q.n_states,q.n_age);
St0(1,:)=0.79*population-10; %Susceptible 1
St0(2,:)=0; %Initially vaccinated
St0(3,:)=10*ones(1,q.n_age); %Infectious 1s
St0(4,:)=0; %Treatment
St0(5,:)=0.2*population; %Recovered
St0(6,:)=0; %Infectious 1r
St0(7,:)=0.01*population; %Carrier sensitive
St0(8,:)=0; %Carrier resistant
St0(9,:)=0; %Susceptible 2
St0(10,:)=0; %Infectious 2s
St0(11,:)=0; %Infectious 2r
%All Initial Vax states are 0

St0 = reshape(St0,q.n_states*q.n_age,1);


[time, St]=ode45(@(time,St)sir_amr_age(time,St,q,v,mass_vac),1:tf,St0);
%%
time(1:t0,:)=[]; %discard the burn-in period for the time variable
St(1:t0,:)=[]; %discard the burn-in period for the State variable
%%
%reshape the vector so age becomes a separate dimension
test = reshape(St, length(St), q.n_states,q.n_age );
St = permute(test,[1 3 2]);
tot_pop_age = sum(St,3);
tot_pop = sum(sum(St,3),2);
S1(:,:,c)=St(:,:,1)+St(:,:,12); %Number fully susceptible through time for control parameter c
S2(:,:,c)=St(:,:,9)+St(:,:,19); %Number susceptible to subclinical infection 
Vacc_age(:,:,c)=St(:,:,2)+St(:,:,22); %Number vaccinated by age
Vacc(:,c)=sum(St(:,:,2),2)+sum(St(:,:,22),2); %Number vaccinated
SensInf(:,:,c)=St(:,:,3)+St(:,:,13); %Number of first sensitive infections (because second are subclinical)
ResInf(:,:,c)=St(:,:,6)+St(:,:,16); %Number of first resistant infections (because second are subclinical)
Inf(:,:,c)=sum(St(:,:,[3,6,10,11,13,16,20,21]),3); %Total number infectious 
Rec(:,:,c)=St(:,:,5)+St(:,:,15); %Number recovered 
SensCar_age(:,:,c)=(St(:,:,7)+St(:,:,17)); %Number of carriers of sensitive strains
ResCar_age(:,:,c)=St(:,:,8)+St(:,:,18); %Number of carriers of resistant strains
Car_age(:,:,c)=St(:,:,7)+St(:,:,8)+St(:,:,17)+St(:,:,18); %Number of total carriers
SensCar(:,c)=sum(St(:,:,7)+St(:,:,17),2); %Number of carriers of sensitive strains
ResCar(:,c)=sum(St(:,:,8)+St(:,:,18),2); %Number of carriers of resistant strains
Car(:,c)=sum(St(:,:,7)+St(:,:,8)+St(:,:,17)+St(:,:,18),2); %Number of total carriers

%assign relative risk vector
rr_vec=ones(1,q.n_age);
rr_vec(1) = q.rr_age1; %.37;
rr_vec(2) = q.rr_age1; %.68;
rr_vec(3) =  q.rr_age2;

lambda_s(:,c)=q.betak*sum(St(:,:,3)+St(:,:,10)+q.rc*St(:,:,7)+St(:,:,13)+St(:,:,20)+q.rc*St(:,:,17),2)/q.N;
lambda_r(:,c)=q.betak*sum(St(:,:,6)+St(:,:,11)+q.rc*St(:,:,8)+St(:,:,16)+St(:,:,21)+q.rc*St(:,:,18),2)/q.N;
Cases_age(:,:,c)=  q.s*(lambda_s(:,c)+lambda_r(:,c)).*repmat(rr_vec,length(St),1).*S1(:,:,c) + q.sv*(1-q.ep)*(lambda_s(:,c)+lambda_r(:,c)).*repmat(rr_vec,length(St),1).*Vacc_age(:,:,c); %number of new cases (who just became infected) at each time step
Cases(:,c)=q.s*(lambda_s(:,c)+lambda_r(:,c)).*sum(S1(:,:,c).*repmat(rr_vec,length(St),1),2) + q.sv*(1-q.ep)*(lambda_s(:,c)+lambda_r(:,c)).*sum(repmat(rr_vec,length(St),1).*Vacc_age(:,:,c),2); %number of new cases (who just became infected) at each time step
Doses_age(:,:,c)=  q.tau*Cases_age(:,:,c); %number of AB Doses (who just got treated) at each time step by age
Doses(:,c)=q.tau*Cases(:,c); %number of Doses at each time step

AnnualInc_age(:,:,c)=Cases_age(:,:,c)*52./(tot_pop_age/1e5); %annual incidence per 100,000 (divide by 10 pop and multiply by 52 weeks)
AnnualInc(:,c)=Cases(:,c)*52./(tot_pop/1e5); %annual incidence per 100,000 (divide by 10 pop and multiply by 52 weeks)
SensCases_age(:,:,c)=q.s*(lambda_s(:,c)).*repmat(rr_vec,length(St),1).*S1(:,:,c) + q.sv*(1-q.ep)*(lambda_s(:,c)).*repmat(rr_vec,length(St),1).*Vacc_age(:,:,c);
ResCases_age(:,:,c)=q.s*(lambda_r(:,c)).*repmat(rr_vec,length(St),1).*S1(:,:,c) + q.sv*(1-q.ep)*(lambda_r(:,c)).*repmat(rr_vec,length(St),1).*Vacc_age(:,:,c);
SensCases(:,c)=q.s*(lambda_s(:,c)).*sum(S1(:,:,c).*repmat(rr_vec,length(St),1),2) + q.sv*(1-q.ep)*(lambda_s(:,c)).*sum(repmat(rr_vec,length(St),1).*Vacc_age(:,:,c),2);
ResCases(:,c)=q.s*(lambda_r(:,c)).*sum(S1(:,:,c).*repmat(rr_vec,length(St),1),2) + q.sv*(1-q.ep)*(lambda_r(:,c)).*sum(repmat(rr_vec,length(St),1).*Vacc_age(:,:,c),2);
AnnualIncSens(:,c)= SensCases(:,c)*52./(tot_pop/1e5);
AnnualIncRes(:,c) = ResCases(:,c)*52./(tot_pop/1e5);
DeathsSens_age(:,:,c) = q.d1s*q.a*(St(:,:,3)+St(:,:,13));
DeathsRes_age(:,:,c)= q.d1r*q.a*(St(:,:,6)+St(:,:,16));
DeathsSens(:,c) = sum(DeathsSens_age(:,:,c),2);
DeathsRes(:,c) = sum(DeathsRes_age(:,:,c),2);
Deaths(:,c) = DeathsSens(:,c) + DeathsRes(:,c);

%Get flow into resistant compartment from de novo vs transmitted:
% Create new matrix of sum of DN and Trans for denominator
ResCasesDN_age(:,:,c)=q.r*St(:,:,4)+q.r*St(:,:,14);
ResCasesDN(:,c)=sum(ResCasesDN_age(:,:,c),2);
ResCasesTrDN_age(:,:,c)= q.s*(lambda_r(:,c)).*repmat(rr_vec,length(St),1).*S1(:,:,c) + q.sv*(1-q.ep)*(lambda_r(:,c)).*repmat(rr_vec,length(St),1).*Vacc_age(:,:,c) +q.r*St(:,:,4)+q.r*St(:,:,14);                       
ResCasesTrDN(:,c)=q.s*(lambda_r(:,c)).*sum(S1(:,:,c).*repmat(rr_vec,length(St),1),2) + q.sv*(1-q.ep)*(lambda_r(:,c)).*sum(repmat(rr_vec,length(St),1).*Vacc_age(:,:,c),2)+ sum(q.r*St(:,:,4)+q.r*St(:,:,14),2);

%%
TotCases=sum(Cases(tvacc+1:end,1));

    Percent_ResCases(:,c)=mean(ResCases(tvacc+1:end,c)./Cases(tvacc+1:end,c));
    Percent_ResCasesDN(:,c)=mean(ResCasesDN(tvacc+1:end,c)./Cases(tvacc+1:end,c));
    Percent_ResCasesTrDN(:,c)=mean(ResCasesTrDN(tvacc+1:end,c)./Cases(tvacc+1:end,c));

    Percent_SensCases(:,c)=mean(SensCases(tvacc+1:end,c)./Cases(tvacc+1:end,c));

    Percent_ResCar(:,c)=ResCar(end,c)/tot_pop(end);
    Percent_SensCar(:,c)=SensCar(end,c)/tot_pop(end);
    TotCasesAverted(:,c)=sum(Cases(tvacc+1:end,1))-sum(Cases(tvacc+1:end,c));
    ResCasesAverted(:,c)=sum(ResCases(tvacc+1:end,1))-sum(ResCases(tvacc+1:end,c));
    Percent_TotCasesAverted(:,c)=(sum(Cases(tvacc+1:end,1))-sum(Cases(tvacc+1:end,c)))./sum(Cases(tvacc+1:end,1));
    Percent_ResCasesAverted(:,c)=(sum(ResCases(tvacc+1:end,1))-sum(ResCases(tvacc+1:end,c)))./sum(ResCases(tvacc+1:end,1));
    ResCar_PrevChange(:,c)= ResCar(end,1)-ResCar(end,c);
    Percent_ResCar_PrevChange(:,c)= (ResCar(end,1)-ResCar(end,c))./ResCar(end,1);
    %by age
    Percent_ResCases_age(:,:,c)=sum(ResCases_age(tvacc+1:end,:,c))./sum(Cases_age(tvacc+1:end,:,c));
    Percent_SensCases_age(:,:,c)=sum(SensCases_age(tvacc+1:end,:,c))./sum(Cases_age(tvacc+1:end,:,c));
    Percent_ResCar_age(:,:,c)=ResCar_age(end,:,c)./tot_pop_age(end,:);
    Percent_SensCar_age(:,:,c)=SensCar_age(end,:,c)/tot_pop(end,:);
    
    DeathsAverted(:,c) = sum(Deaths(tvacc+1:end,1))-sum(Deaths(tvacc+1:end,c));
    Percent_DeathsAverted(:,c)=( sum(Deaths(tvacc+1:end,1))-sum(Deaths(tvacc+1:end,c)))./sum(Deaths(tvacc+1:end,1));
    
    DosesAverted(:,c) = sum(Doses(tvacc+1:end,1))-sum(Doses(tvacc+1:end,c));
    Percent_DosesAverted(:,c) = (sum(Doses(tvacc+1:end,1))-sum(Doses(tvacc+1:end,c)))./sum(Doses(tvacc+1:end,1));
    
    PercentResCasesMat(j,k,l,m,c) = Percent_ResCases(:,c);
    PercentResCasesMatDN(j,k,l,m,c) = Percent_ResCasesDN(:,c);
    PercentResCasesMatTrDN(j,k,l,m,c) = Percent_ResCasesTrDN(:,c);
    PercentSensCasesMat(j,k,l,m,c) = Percent_SensCases(:,c);
    

                end
            end
        end
    end
end
save('R0Outputs')
%% Calculate difference in percentage resistant between no vacc and camp15
VaccStrats = {'No Vacc','Routine','Camp 5','Camp 15'};
DiffResC15_pc = (PercentResCasesMat(:,:,:,:,1)-PercentResCasesMat(:,:,:,:,4))./PercentResCasesMat(:,:,:,:,1);
DiffResC15 = (PercentResCasesMat(:,:,:,:,1)-PercentResCasesMat(:,:,:,:,4));
%% Do Scatter plots of change in percentage resistant vs R0

% Create vector of Diffs unwrapped by R0
R0DiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:length(R0samp)
    %get indices to read in
    tmpIdx = [(i-1)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2)+1:i*size(sympSamp,2)*size(CovMat,1)*size(resParams,2)];
    %read in R0 value
    R0DiffMat(tmpIdx,1)=R0samp(i);
    %read in Diff value
    R0DiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(i,:,:,:),size(sympSamp,2)*size(CovMat,1)*size(resParams,2),1);
end
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,1)
scatter(R0DiffMat(:,1),R0DiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("R_0 Value")
ytickformat(gca, 'percentage');
title("A)")

[corrR0,pR0] =corr(R0DiffMat(:,1),R0DiffMat(:,2));
%% Do Scatter plots of change in percentage resistant vs prop symptomatic

% Create vector of Diffs unwrapped by prop symptomatic
sympDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:length(sympSamp)
    %get indices to read in
    tmpIdx = [(i-1)*size(R0samp,2)*size(CovMat,1)*size(resParams,2)+1:i*size(R0samp,2)*size(CovMat,1)*size(resParams,2)];
    %read in symp value
    sympDiffMat(tmpIdx,1)=sympSamp(i);
    %read in Diff value
    sympDiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(:,i,:,:),size(R0samp,2)*size(CovMat,1)*size(resParams,2),1);
end
subplot(2,4,2)
scatter(sympDiffMat(:,1),sympDiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("Proportion Symptomatic")
ytickformat(gca, 'percentage');
title("B)")
[corrSymp,pSymp] =corr(sympDiffMat(:,1),sympDiffMat(:,2));
%% Do Scatter plots of change in percentage resistant vs coverage

% Create vector of Diffs unwrapped by Difference in Coverage
covDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:size(CovMat,1)
    %get indices to read in
    tmpIdx = [(i-1)*size(R0samp,2)*size(sympSamp,2)*size(resParams,2)+1:i*size(R0samp,2)*size(sympSamp,2)*size(resParams,2)];
    %read in R0 value
    covDiffMat(tmpIdx,1)=100*mean(CovMat(i,:));
    %read in Diff value
    covDiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(:,:,i,:),size(R0samp,2)*size(sympSamp,2)*size(resParams,2),1);
end
subplot(2,4,3)
scatter(covDiffMat(:,1),covDiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("Mean Routine Vaccine Coverage")
ytickformat(gca, 'percentage');
xtickformat(gca, 'percentage');
set(gca,'XLim',[50 100])
title("C)")        
[corrCov,pCov] =corr(covDiffMat(:,1),covDiffMat(:,2));
%% Test change in proportion resistant vs resistance parameter values

%% vs r, Rate of resistance acquisition from treatment
% Create vector of Diffs unwrapped by resParams
rDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:size(resParams,2)
    %get indices to read in
    tmpIdx = [(i-1)*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)+1:i*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)];
    %read in R0 value
    rDiffMat(tmpIdx,1)=resParams(i).r;
    %read in Diff value
    rDiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(:,:,:,i),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1),1);
end
subplot(2,4,5)
scatter(rDiffMat(:,1),rDiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("Rate of resistance acquisition from treatment")
ytickformat(gca, 'percentage');
title("E)")
        
[corrr,pr] =corr(rDiffMat(:,1),rDiffMat(:,2));

%% vs d1r, rate of recovery from primary infection with a resistant strain
d1rDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:size(resParams,2)
    %get indices to read in
    tmpIdx = [(i-1)*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)+1:i*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)];
    %read in R0 value
    d1rDiffMat(tmpIdx,1)=resParams(i).d1r;
    %read in Diff value
    d1rDiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(:,:,:,i),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1),1);
end
subplot(2,4,6)
scatter(d1rDiffMat(:,1),d1rDiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("Rate of recovery")
ytickformat(gca, 'percentage');
title("F)")        
[corrd1r,pd1r] =corr(d1rDiffMat(:,1),d1rDiffMat(:,2));

%% vs d2r, rate of recovery from subclinical infection with a resistant strain
d2rDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:size(resParams,2)
    %get indices to read in
    tmpIdx = [(i-1)*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)+1:i*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)];
    %read in R0 value
    d2rDiffMat(tmpIdx,1)=resParams(i).d2r;
    %read in Diff value
    d2rDiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(:,:,:,i),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1),1);
end
subplot(2,4,7)
scatter(d2rDiffMat(:,1),d2rDiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("Rate of recovery")
ytickformat(gca, 'percentage');
title("G)")
%xtickformat(gca, 'percentage');
%set(gca,'XLim',[50 100])
%title(["c) Percent Reduction in Resistance Proportion vs", "Rate of Recovery from Subclinical Resistant Infection"])
%text(.4,55,['corr = ',num2str(corr(d2rDiffMat(:,1),d2rDiffMat(:,2)))])
        
[corrd2r,pd2r] =corr(d2rDiffMat(:,1),d2rDiffMat(:,2));
%% vs rr, Relative risk of transmission for resistant strain
rrDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),2);
for i = 1:size(resParams,2)
    %get indices to read in
    tmpIdx = [(i-1)*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)+1:i*size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)];
    %read in R0 value
    rrDiffMat(tmpIdx,1)=resParams(i).rr;
    %read in Diff value
    rrDiffMat(tmpIdx,2)=100*reshape(DiffResC15_pc(:,:,:,i),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1),1);
end
subplot(2,4,8)
scatter(rrDiffMat(:,1),rrDiffMat(:,2))
ylabel("Percent Reduction in Resistance")
xlabel("Relative Risk")
ytickformat(gca, 'percentage');
title("H)")        
[corrrr,prr] =corr(rrDiffMat(:,1),rrDiffMat(:,2));
%% Get Linear model of percent change in resistance by R0 and symp

% Create vector of Diffs unwrapped by R0 and sympSamp
R0SympDiffMat = zeros(length(R0samp)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),3);
for i = 1:length(R0samp)
    %get indices to read in R0
    tmpIdx = [(i-1)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2)+1:i*size(sympSamp,2)*size(CovMat,1)*size(resParams,2)];
    %read in R0 value
    R0SympDiffMat(tmpIdx,1)=R0samp(i);
    for j = 1:size(sympSamp,2)
        %get indices to read in symp
        tmpIdx2 = (tmpIdx(1)-1)+[(j-1)*size(CovMat,1)*size(resParams,2)+1:j*size(CovMat,1)*size(resParams,2)];

         %read in symp value
         R0SympDiffMat(tmpIdx2,2) = sympSamp(j);
         %read in Diff value
         R0SympDiffMat(tmpIdx2,3)=100*reshape(DiffResC15_pc(i,j,:,:),size(CovMat,1)*size(resParams,2),1);

    end
end
%Create table
tblR0SympDiff= array2table(R0SympDiffMat,'VariableNames',{'R0','Symp','ResDiff'});
resDiff_mdl = fitlm(tblR0SympDiff,'ResDiff~R0+Symp');

tblR0SympDiff.InitRes=reshape(PercentResCasesMat(:,:,:,:,1),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),1);
resDiff_mdl_1 = fitlm(tblR0SympDiff,'ResDiff~R0+Symp+InitRes');


%% Get range of changes in Resistant proportion for different starting resistance percentages
ResPCVec = reshape(PercentResCasesMat(:,:,:,:,1),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),1);
DiffVec = reshape(-DiffResC15_pc,size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),1);
%cateogorize resistance percentages
% define bins for starting resistance percentages
resBins= [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1];
ResCat = discretize(ResPCVec,resBins,resBins(1:end-1));
%figure
subplot(2,4,4)
boxplot(-100*DiffVec,ResCat)
ytickformat(gca,'percentage')
set(gca,'FontSize',12,'YLim',[0 60])
xlabel("Initial Proportion Resistant",'FontSize',13)
title("D)")
saveas(gcf,[Folder,'/FigS7_ResDiffCorrFigs'],'jpg');

%% Output
[DiffMedians, DiffCIs]=grpstats(DiffVec,ResCat,{'median','predci'},'Alpha',.05);
DiffTable=table(resBins(1:end-1)',DiffMedians, DiffCIs(:,1),DiffCIs(:,2),'VariableNames',{'InitialRes','MedianDiff','LowCI','HighCI'});
save('SampleResPCs','DiffTable','ResPCVec','DiffVec','resParams','ResCat','tblR0SympDiff','resDiff_mdl');
%% Plot difference in proportion of cases from transmitted res vs de novo res in BL vs RC15
%divide flow into res from transmitted cases by total flow (trans plus de
%novo)
propTrans=PercentResCasesMat./PercentResCasesMatTrDN;
% get proportions for BL and RC15
propTransBL = reshape(propTrans(:,:,:,:,1),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),1);
propTransRC15 = reshape(propTrans(:,:,:,:,4),size(R0samp,2)*size(sympSamp,2)*size(CovMat,1)*size(resParams,2),1);

figure
boxplot(100*[propTransBL,propTransRC15], 'Colors','k','Labels',{'Baseline','Camp15'})
ytickformat(gca, 'percentage');
title('Percentage of Resistant Cases  Over 10 Years from Transmitted Resistance')
saveas(gcf,[Folder,'/FigS9_PercCasesTrans'],'jpg');
