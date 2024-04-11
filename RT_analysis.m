%% busca aca C:\Users\ALEJANDRO\Desktop\FIGURES_PAPER1\Figs\Fig2\conAMA\conAMA070523
clear all
load RT_conflict_effect.mat
% trials ordered as:
% incongruent awake, congruent awake, incongruent drowsy, congruent drowsy.
Nsub = length(RT_conflict_effect{1});

%% Load Data
sub_dat = cell(1,Nsub);
for subi=1:Nsub
    inc_aw = RT_conflict_effect{1}{subi}(:);
    con_aw = RT_conflict_effect{2}{subi}(:);
    inc_dr = RT_conflict_effect{3}{subi}(:);
    con_dr = RT_conflict_effect{4}{subi}(:);
    RT = [inc_aw; con_aw; inc_dr; con_dr];
    % congruence coded as 0: congruent, 1: incongruent
    congruent = [ones(size(inc_aw)); zeros(size(con_aw)); ones(size(inc_dr)); zeros(size(con_dr))];
    % drowsy coded as 0: alert, 1: drowsy
    drowsy = [zeros(size(inc_aw)); zeros(size(con_aw)); ones(size(inc_dr)); ones(size(con_dr))];
    dat = table(RT,congruent,drowsy);
    dat.ID = repmat(subi,length(RT),1);
    sub_dat{subi} = dat;
end
% concatenate all data into a table
RT_dat = vertcat(sub_dat{:});
RT_dat.logRT = log(RT_dat.RT);

%
%
%

%% fit group mixed-effect model
mdl = fitlme(RT_dat,'RT~congruent*drowsy + (congruent * drowsy|ID)');
disp(mdl)

%
%
%

%% fit within participant models
sub_mdls = cell(1,Nsub);
for subi=1:Nsub
    sub_dat = RT_dat(RT_dat.ID==subi,:);
    mdl = fitlm(sub_dat, 'RT ~ congruent*drowsy');
    sub_mdls{subi} = mdl;
end
pvals = cell2mat(cellfun(@(x) x.Coefficients.pValue, sub_mdls,'UniformOutput',false));
estimates = cell2mat(cellfun(@(x) x.Coefficients.Estimate, sub_mdls,'UniformOutput',false));

% coefficients stored in axis 1 as:
% 1: intercept, 2: congruence, 3:, drowsy, 4: interaction
fprintf(1,"\nNumber of significant congruence coefficients: %d\n", sum(pvals(2,:)<0.05));
fprintf(1,"Number of significant drowsy coefficients: %d\n", sum(pvals(3,:)<0.05));
fprintf(1,"Number of significant interaction coefficients: %d\n", sum(pvals(4,:)<0.05));

%
map = bayesprev_map(sum(pvals(2,:)<0.05), Nsub, 0.05)
%needs optimization toolbox
hpdi = bayesprev_hpdi(0.95, sum(pvals(2,:)<0.05), Nsub, 0.05)

%
% example of Bayesian prevalence analyses on a single plot for congruency
figure
k = sum(pvals(2,:)<0.05);
n = Nsub;
co = get(gca,'ColorOrder'); ci=1;
hold on

x = linspace(0,1,100);
posterior = bayesprev_posterior(x,k,n,alpha);
plot(x, posterior,'Color',co(ci,:));

% add MAP as a point
xmap = bayesprev_map(k,n,alpha);
pmap = bayesprev_posterior(xmap,k,n,alpha);
plot(xmap, pmap,'.','MarkerSize',20,'Color',co(ci,:));

% add lower bound as a vertical line
bound = bayesprev_bound(0.95,k,n,alpha);
line([bound bound], [0 bayesprev_posterior(bound,k,n,alpha)],'Color',co(ci,:),'LineStyle',':')

% add 95% HPDI
oil = 2;
iil = 4;
h = bayesprev_hpdi(0.95,k,n,alpha);
plot([h(1) h(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',oil)
% add 50% HPDI
h = bayesprev_hpdi(0.5,k,n,alpha);
plot([h(1) h(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',iil)

xlabel('Population prevalence proportion')
ylabel('Posterior density')

%

%%  MI analysis 
Nbin = 3;
Irt_cong = zeros(1,Nsub);
Irt_drowsy = zeros(1,Nsub);
coI_cong_drowsy = zeros(1,Nsub);
Irt_cong_sig = zeros(1,Nsub);
Irt_drowsy_sig = zeros(1,Nsub);
coI_cong_drowsy_sig = zeros(1,Nsub);
Nperm = 1000;
Irt_cong_perm = zeros(1,Nperm);
Irt_drowsy_perm = zeros(1,Nperm);
coI_cong_drowsy_perm = zeros(1,Nperm);

Irt_cong_aw = zeros(1,Nsub);
Irt_cong_dr = zeros(1,Nsub);
for subi=1:Nsub
    subidx = RT_dat.ID == subi;
    sub_dat = RT_dat(subidx,:);
    congruent = sub_dat.congruent;
    drowsy = sub_dat.drowsy;
    RT = sub_dat.RT;
    % bin each alertness condition separately
    bRT = [eqpopbin(RT(drowsy==0), Nbin); eqpopbin(RT(drowsy==1),Nbin)];
    % bin all RT's together
%     bRT = eqpopbin(RT, Nbin)
    Irt_cong(subi) = calcinfo(bRT, Nbin, congruent, 2);
    Irt_drowsy(subi) = calcinfo(bRT, Nbin, drowsy, 2);
    coI_cong_drowsy(subi) = Irt_cong_sig(subi) + Irt_drowsy(subi) - ...
                            calcinfo(bRT,Nbin,drowsy*2+congruent,4);

    for pi=1:Nperm
        idx = randperm(length(bRT));
        Irt_cong_perm(pi) = calcinfo(bRT(idx),Nbin,congruent,2);
        Irt_drowsy_perm(pi) = calcinfo(bRT,Nbin,drowsy(idx),2);

        % permute for coI. keep relationship with congruence, shuffle
        % drowsy condition
        coI_cong_drowsy_perm(pi) = Irt_cong(subi) + Irt_drowsy_perm(pi) - ...
                            calcinfo(bRT, Nbin, drowsy(idx)*2+congruent, 4);
    end
    Irt_cong_sig(subi) = Irt_cong(subi)>prctile(Irt_cong_perm,95);
    Irt_drowsy_sig(subi) = Irt_drowsy(subi)>prctile(Irt_drowsy_perm,95);
    % one sided test for negative coI 
    coI_cong_drowsy_sig(subi) = coI_cong_drowsy(subi)<prctile(coI_cong_drowsy_perm,5);

    idx = drowsy==0;
    bRT_cond = eqpopbin(RT(idx),Nbin);
    Irt_cong_aw(subi) = calcinfo(bRT_cond, Nbin, congruent(idx), 2);
    idx = drowsy==1;
    bRT_cond = eqpopbin(RT(idx),Nbin);
    Irt_cong_dr(subi) = calcinfo(bRT_cond, Nbin, congruent(idx), 2);
end
% per-subject difference in congruence-RT effect between awake and drowsy
I_drowsy_diff = Irt_cong_aw - Irt_cong_dr;
%%
fprintf(1,"\nNumber of significant I(RT; cong): %d\n", sum(Irt_cong_sig));
fprintf(1,"Number of significant drowsy coefficients: %d\n", sum(Irt_drowsy_sig));
fprintf(1,"Number of significant interaction coefficients: %d\n", sum(coI_cong_drowsy_sig));

%% 
% Estimate prevalence at
% https://estimate.prevalence.online/
% or use code from
% https://github.com/robince/bayesian-prevalence

map = bayesprev_map(sum(coI_cong_drowsy_sig), Nsub, 0.05)
%needs optimization toolbox
hpdi = bayesprev_hpdi(0.95, sum(coI_cong_drowsy_sig), Nsub, 0.05)
