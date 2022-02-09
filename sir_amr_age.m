function dSt=sir_amr_age(t,St_vec,q,v,mass_vac) %,fitted,unkp)

if nargin < 6
    fitted = 0;
    unkp = 0;
end
%set up passed in parameters
n_age = q.n_age;
n_states = q.n_states;
mub = q.mub;
mu = q.mu;
d1s = q.d1s;
a = q.a;
w2 = q.w2;
t2s = q.t2s;
t2r = q.t2r;
betak = q.betak;
e = q.e;
e_i = q.e_i;
ep = q.ep;
w1 = q.w1;
s = q.s;
sv = q.sv;
tau = q.tau;
g = q.g;
r = q.r;
d1r = q.d1r;
rr = q.rr;
rc = q.rc;
t1s_y = q.t1s_y;
t1s_o = q.t1s_o;
t1r_y = q.t1r_y;
t1r_o = q.t1r_o;
rd = q.rd;
d2s = q.d2s;
d2r = q.d2r;
rr_age1 = q.rr_age1;
rr_age2 = q.rr_age2;

%try setting resistance to not emerge until 10 years before vaccination
%intervention
if t <4680 % 5 years =4940%20 years:4160 %10 years before intervention: 4680
    r=0;
end

% % check if fitting is happening and assign new values if so
% if fitted ==1
%     betak = unkp.betak; %effective contact rate
%     rr_age1 = unkp.rr_age1; %relative risk for age group 1
%     rr_age2 = unkp.rr_age2; %rel risk for age group 2
%     r = unkp.r; %%treatment-induced/acquired resistance rate
%     s = unkp.s; %proportion symptomatic
% end
    


%Reshape for 6 age groups: 0-9 months, 9 months - 2 years, 2-5, 5-15, 15-25, 25- end
St= reshape(St_vec,n_states,n_age);

%Defining state variables
S1=St(1,:);
V1=St(2,:);
I1s=St(3,:);
It=St(4,:);
R=St(5,:);
I1r=St(6,:);
Cs=St(7,:);
Cr=St(8,:);
S2=St(9,:);
I2s=St(10,:);
I2r=St(11,:);
% for counting vaccinations
S1_v = St(12,:);
I1s_v = St(13,:);
It_v = St(14,:);
R_v = St(15,:);
I1r_v = St(16,:);
Cs_v = St(17,:);
Cr_v = St(18,:);
S2_v = St(19,:);
I2s_v=St(20,:);
I2r_v=St(21,:);
% Second vaccination compartment
V2 = St(22,:);
N = sum(sum(St));

% look into this
%betap = R0*(delta+mean(mu))/(1 + rC*delta*(agep*(theta/mean(mu)')))*[b1*ones(1,al); b2*ones(1,al); ones(al-2,al)]*al/max(eig([b1*ones(1,al); b2*ones(1,al); ones(al-2,al)])); %person-to-person transmission rate
%lambdap=betap*(y(al+1:2*al)+r*y(4*al+1:5*al)+rC*y(5*al+1:6*al))/sum(y(1:6*al));


lambdas=betak*(sum(I1s)+sum(I2s)+rc*sum(Cs)+sum(I1s_v)+sum(I2s_v)+rc*sum(Cs_v))/N;
lambdar=betak*rr*(sum(I1r)+sum(I2r)+rc*sum(Cr)+sum(I1r_v)+sum(I2r_v)+rc*sum(Cr_v))/N;

%assign relative risk vector
rr_vec=ones(1,n_age);
rr_vec(1) = rr_age1; %.37;
rr_vec(2) = rr_age1; %rr_age2 %.68;
rr_vec(3) = rr_age2;

%assign theta vectors
t1s = t1s_y*ones(1,n_age); %young proportion carriers
t1r = t1r_y*ones(1,n_age);
t1s(n_age) = t1s_o; %old proportion carriers
t1r(n_age) = t1r_o;

%set up aging vectors
%age_rate = [1/(0.75*52); 1/(1.25*52); 1/(3*52); 1/520; 1/520; 0];
%above from original age groups
age_rate = q.age_rate; %[1/(0.75*52); 1/(1.25*52); 1/(3*52); 1/260; 1/260; 1/780; 0];
%aging in
u_in = zeros(n_states,n_age);
u_in(1,1)= mub*N +a*d1s*sum(I1s+I1s_v)+a*d1r*sum(I1r+I1r_v);
%matrix of aging rates times disease states from age class 1 younger
u_in(:,2:n_age) =St(:,1:n_age-1).*repmat(age_rate(1:n_age-1)',n_states,1);

%aging out
%aging rates times current age states
u_out = St.*repmat(age_rate', n_states,1); 

% Set up vaccination vectors
% for routine: only people in first age group get vaccinated as they age
% into second age group
vacc_rout = zeros(1,n_age);
vacc_rout(2) = v(round(t));

% for 5 yr catch-up campaign: people in second and third age groups get vaccinated for 4 weeks
% for 15 year catch-up campaign: people in second, third and fourth age
% groups get vaccinated for 4 weeks

vacc_camp = mass_vac(round(t),:);

%Differential equations 
% Basic equations
% dSt(1,:)=u_in(1,:) - u_out(1,:) - (lambdas+lambdar)*S1 - mu*S1; %dS1/dt
% dSt(2,:)=u_in(2,:) - u_out(2,:) - w1*V - (1-ep)*(lambdas+lambdar)*V - mu*V; %dV/dt
% dSt(3,:)=u_in(3,:) - u_out(3,:) +lambdas*(1-s*tau)*S1 - d1s*I1s - mu*I1s; %dI1s/dt
% dSt(4,:)=u_in(4,:) - u_out(4,:) + lambdas*s*tau*S1 - g*It - r*It - mu*It; %dIt/dt
% dSt(5,:)=u_in(5,:) - u_out(5,:) + d1s*(1-a-t1s)*I1s + d1r*(1-a-t1r)*I1r + d2s*(1-t2s)*I2s + d2r*(1-t2r)*I2r + g*It - w2*R - mu*R; %dR/dt
% dSt(6,:)=u_in(6,:) - u_out(6,:) + lambdar*S1 - d1r*I1r + r*It - mu*I1r; %dI1r/dt
% dSt(7,:)=u_in(7,:) - u_out(7,:) + d1s*t1s*I1s + d2s*t2s*I2s - mu*Cs; %dCs/dt
% dSt(8,:)=u_in(8,:) - u_out(8,:) + d1r*t1r*I1r + d2r*t2r*I2r - mu*Cr; %dCr/dt
% dSt(9,:)=u_in(9,:) - u_out(9,:) + w2*R - (lambdas+lambdar)*S2 - mu*S2; %dS2/dt
% dSt(10,:)=u_in(10,:) - u_out(10,:) + lambdas*S2 - d2s*I2s - mu*I2s; %dI2s/dt
% dSt(11,:)=u_in(11,:) - u_out(11,:) + lambdar*S2 - d2r*I2r - mu*I2r; %dI2r/dt
% % Vaccinated infected individuals (to keep count of vaccinations)
% % disease process for individuals who have been vaccinated before.
% dSt(12,:) = u_in(12,:) - u_out(12,:) + w1*V - (lambdas+lambdar)*S1_v - mu*S1_v; %dS1_v/dt -- waning immunity from vaccination goes here
% dSt(13,:)=u_in(13,:) - u_out(13,:) +lambdas*((1-s*tau)*S1_v + (1-ep)*(1-sv*tau)*V)- d1s*I1s_v - mu*I1s_v; %dI1s_v/dt people infected bc leaky vaccine go here
% dSt(14,:)=u_in(14,:) - u_out(14,:) + lambdas*(s*tau*S1_v + (1-ep)*sv*tau*V) - g*It_v - r*It_v - mu*It_v; %dIt/dt
% dSt(15,:)=u_in(15,:) - u_out(15,:) + d1s*(1-a-t1s)*I1s_v + d1r*(1-a-t1r)*I1r_v + d2s*(1-t2s)*I2s_v + d2r*(1-t2r)*I2r_v + g*It_v - w2*R_v - mu*R_v; %dR/dt
% dSt(16,:)=u_in(16,:) - u_out(16,:) + lambdar*(S1_v+(1-ep)*V) - d1r*I1r_v + r*It_v - mu*I1r_v; %dI1r/dt
% dSt(17,:)=u_in(17,:) - u_out(17,:) + d1s*t1s*I1s_v + d2s*t2s*I2s_v - mu*Cs_v; %dCs/dt
% dSt(18,:)=u_in(18,:) - u_out(18,:) + d1r*t1r*I1r_v + d2r*t2r*I2r_v - mu*Cr_v; %dCr/dt
% dSt(19,:)=u_in(19,:) - u_out(19,:) + w2*R_v - (lambdas+lambdar)*S2_v - mu*S2_v; %dS2/dt
% dSt(20,:)=u_in(20,:) - u_out(20,:) + lambdas*S2_v - d2s*I2s_v - mu*I2s_v; %dI2s/dt
% dSt(21,:)=u_in(21,:) - u_out(21,:) + lambdar*S2_v - d2r*I2r_v - mu*I2r_v; %dI2r/dt

%Add in vaccination
%disp(mu);

dSt(1,:)=u_in(1,:) - u_out(1,:) - (lambdas+lambdar)*rr_vec.*S1 - mu.*S1; %dS1/dt
dSt(2,:)=u_in(2,:) - u_out(2,:) - w1*V1 - (1-ep)*(lambdas+lambdar)*rr_vec.*V1 - mu.*V1; %dV1/dt
dSt(3,:)=u_in(3,:) - u_out(3,:) +lambdas*(1-s*tau)*rr_vec.*S1 - d1s*I1s - mu.*I1s; %dI1s/dt
dSt(4,:)=u_in(4,:) - u_out(4,:) + lambdas*s*tau*rr_vec.*S1 - g*It - r*It - mu.*It; %dIt/dt
dSt(5,:)=u_in(5,:) - u_out(5,:) + d1s*(-t1s+1-a).*I1s + d1r*(-t1r+1-a).*I1r + d2s*(1-t2s)*I2s + d2r*(1-t2r)*I2r + g*It - w2*R - mu.*R; %dR/dt
dSt(6,:)=u_in(6,:) - u_out(6,:) + lambdar*rr_vec.*S1 - d1r*I1r + r*It - mu.*I1r; %dI1r/dt
dSt(7,:)=u_in(7,:) - u_out(7,:) + d1s*t1s.*I1s + d2s*t2s*I2s - mu.*Cs; %dCs/dt
dSt(8,:)=u_in(8,:) - u_out(8,:) + d1r*t1r.*I1r + d2r*t2r*I2r - mu.*Cr; %dCr/dt
dSt(9,:)=u_in(9,:) - u_out(9,:) + w2*R +w1*V2- (lambdas+lambdar)*rr_vec.*S2 - mu.*S2; %dS2/dt
dSt(10,:)=u_in(10,:) - u_out(10,:) + lambdas*rr_vec.*S2 - d2s*I2s - mu.*I2s; %dI2s/dt
dSt(11,:)=u_in(11,:) - u_out(11,:) + lambdar*rr_vec.*S2 - d2r*I2r - mu.*I2r; %dI2r/dt
% Vaccinated infected individuals (to keep count of vaccinations)
% disease process for individuals who have been vaccinated before.
dSt(12,:) = u_in(12,:) - u_out(12,:) + w1*V1 - (lambdas+lambdar)*rr_vec.*S1_v - mu.*S1_v; %dS1_v/dt -- waning immunity from vaccination goes here
dSt(13,:)=u_in(13,:) - u_out(13,:) +lambdas*((1-s*tau)*rr_vec.*S1_v + (1-ep)*(1-sv*tau)*rr_vec.*V1)- d1s*I1s_v - mu.*I1s_v; %dI1s_v/dt people infected bc leaky vaccine go here
dSt(14,:)=u_in(14,:) - u_out(14,:) + lambdas*rr_vec.*(s*tau*S1_v + (1-ep)*sv*tau*V1) - g*It_v - r*It_v - mu.*It_v; %dIt/dt
dSt(15,:)=u_in(15,:) - u_out(15,:) + d1s*(-t1s+1-a).*I1s_v + d1r*(-t1r+1-a).*I1r_v + d2s*(1-t2s)*I2s_v + d2r*(1-t2r)*I2r_v + g*It_v - w2*R_v - mu.*R_v; %dR/dt
dSt(16,:)=u_in(16,:) - u_out(16,:) + lambdar*rr_vec.*(S1_v+(1-ep)*V1) - d1r*I1r_v + r*It_v - mu.*I1r_v; %dI1r/dt
dSt(17,:)=u_in(17,:) - u_out(17,:) + d1s*t1s.*I1s_v + d2s*t2s*I2s_v - mu.*Cs_v; %dCs/dt
dSt(18,:)=u_in(18,:) - u_out(18,:) + d1r*t1r.*I1r_v + d2r*t2r*I2r_v - mu.*Cr_v; %dCr/dt
dSt(19,:)=u_in(19,:) - u_out(19,:) + w2*R_v - (lambdas+lambdar)*rr_vec.*S2_v - mu.*S2_v; %dS2/dt
dSt(20,:)=u_in(20,:) - u_out(20,:) + lambdas*rr_vec.*S2_v - d2s*I2s_v - mu.*I2s_v; %dI2s/dt
dSt(21,:)=u_in(21,:) - u_out(21,:) + lambdar*rr_vec.*S2_v - d2r*I2r_v - mu.*I2r_v; %dI2r/dt
% added V2 compartment
dSt(22,:) = u_in(22,:) -u_out(22,:) - w1*V2 - (1-ep)*(lambdas+lambdar)*rr_vec.*V2 - mu.*V2; %dV/dt
%add in routine vaccination if it's time
if v(round(t)) >0
    %take vaccinated out of susceptible and infecteds (already aging out, don't need to subtract again)
    %and add to vaccinated/vaccinated infected
    %Susceptible: go into vaccinated and Susceptible -- take out people who
    %now flow into vaccinated and susceptible vaccinated
    dSt(1,:) = dSt(1,:) - vacc_rout.*u_in(1,:);
    dSt(2,:) = dSt(2,:) + e*vacc_rout.*u_in(1,:);
    dSt(12,:) = dSt(12,:) + (1-e)*vacc_rout.*u_in(1,:);
    
    %For V2, take people out of S2 and R and put them into V2
    dSt(5,:) = dSt(5,:) -vacc_rout.*u_in(5,:);
    dSt(9,:) = dSt(9,:) - vacc_rout.*u_in(9,:);
    dSt(22,:) = dSt(22,:) + e*vacc_rout.*(u_in(5,:)+u_in(9,:));
    dSt(15,:) = dSt(15,:) +(1-e)*vacc_rout.*u_in(15,:);
    dSt(19,:) = dSt(19,:) +(1-e)*vacc_rout.*u_in(19,:);
    
    
    % Infected categories - put into V2
    dSt(3:11,:) = dSt(3:11,:) - vacc_rout.*u_in(3:11,:);
    dSt(22,:) = dSt(22,:) + e_i*sum(vacc_rout.*u_in(3:11,:),1);
    dSt(13:21,:) = dSt(13:21,:) + (1-e_i)*vacc_rout.*u_in(3:11,:);
    
end
   %add in vaccination campaign  
if any(mass_vac(round(t),:))
    %move people into vaccinated status from corresponding age-group
    %Susceptible to vaccinated
    dSt(1,:) = dSt(1,:) - vacc_camp.*S1;
    dSt(2,:) = dSt(2,:) + e*vacc_camp.*S1;
    dSt(12,:) = dSt(12,:) + (1-e)*vacc_camp.*S1;
    %For V2, take people out of S2 and R and put them into V2
    dSt(5,:) = dSt(5,:) -vacc_camp.*u_in(5,:);
    dSt(9,:) = dSt(9,:) - vacc_camp.*u_in(9,:);
    dSt(22,:) = dSt(22,:) + e*vacc_camp.*(u_in(5,:)+u_in(9,:));
    dSt(15,:) = dSt(15,:) +(1-e)*vacc_camp.*u_in(15,:);
    dSt(19,:) = dSt(19,:) +(1-e)*vacc_camp.*u_in(19,:);
    
    % Infected categories
    dSt(3:11,:) = dSt(3:11,:) - vacc_camp.*St(3:11,:);
    dSt(22,:) = dSt(22,:) + e_i*sum(vacc_camp.*St(3:11,:),1);
    dSt(13:21,:) = dSt(13:21,:) + (1-e_i)*vacc_camp.*St(3:11,:);

end
    
 
dSt = reshape(dSt,n_states*n_age,1);
end

