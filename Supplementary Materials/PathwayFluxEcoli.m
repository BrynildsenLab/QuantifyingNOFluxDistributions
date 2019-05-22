function PathwayFluxEcoli(t,dV,model)

EXautox_list = {'EX_NOautox';'EX_NO_NO2r';'EX_N2O3diss'};
EXgas_list = {'EX_NOout'};
autox_list = {'NOautox';'NO_NO2r';'N2O3diss'};
hmp_list = {'HMP6';'HMP11';'HMP13';'HMP14';'HMP16';'HMP17';
       'HMP21';'HMP22';'HMP28';'hmpFe3h_onoo_deg';
       'hmpFe2_no_deg';'hmpFe3_onoo_deg';'hmpFe2h_no_deg';
       'hmpFe3h2_onoo_deg';'hmpFe2h2_no_deg'};
norv_list = {'NORVno'};
fes_list = {'2Fe2S_NO';'4Fe4S_NO'};
nrfa_list = {'NRFAno'};
o2s_list = {'NO_O2S'};




rxnNames = model.rxnNames;
[~,EXautox_ind] = ismember(EXautox_list,rxnNames);
[~,EXgas_ind] = ismember(EXgas_list,rxnNames);
[~,autox_ind] = ismember(autox_list,rxnNames);
[~,hmp_ind] = ismember(hmp_list,rxnNames);
[~,norv_ind] = ismember(norv_list,rxnNames);
[~,fes_ind] = ismember(fes_list,rxnNames);
[~,nrfa_ind] = ismember(nrfa_list,rxnNames);
[~,o2s_ind] = ismember(o2s_list,rxnNames);


cell_frac=model.Frxn(1,1);
media_frac=1-cell_frac;

S=model.S;
[~,NOind] = ismember('no',model.species);

% time-course fluxes through each consumption pathway
V_EXautox = sum(dV(:,EXautox_ind).*repmat(S(NOind,EXautox_ind),size(dV,1),1),2).*media_frac;
V_EXgas = sum(dV(:,EXgas_ind).*repmat(S(NOind,EXgas_ind),size(dV,1),1),2).*media_frac;
V_autox = sum(dV(:,autox_ind).*repmat(S(NOind,autox_ind),size(dV,1),1),2).*cell_frac;
V_hmp = sum(dV(:,hmp_ind).*repmat(S(NOind,hmp_ind),size(dV,1),1),2).*cell_frac;
V_norv = sum(dV(:,norv_ind).*repmat(S(NOind,norv_ind),size(dV,1),1),2).*cell_frac;
V_nrfa = sum(dV(:,nrfa_ind).*repmat(S(NOind,nrfa_ind),size(dV,1),1),2).*cell_frac;
V_fes = sum(dV(:,fes_ind).*repmat(S(NOind,fes_ind),size(dV,1),1),2).*cell_frac;
V_o2s = sum(dV(:,o2s_ind).*repmat(S(NOind,o2s_ind),size(dV,1),1),2).*cell_frac;



Vmat = [V_EXautox,V_EXgas,V_hmp,V_norv,V_nrfa,V_fes,V_o2s,V_autox];
Vmat=abs(Vmat);
Vmat= cumtrapz(t,Vmat);
cellular=Vmat(:,3:end);
Vmat(:,3:end)=[];
Vmat(:,3)=sum(cellular,2);


%Total Cumulative NO Flux Plot
A=area(t,Vmat);
A(1).FaceColor=[0 0 1];
A(2).FaceColor=[1 0 0];
A(3).FaceColor=[1 0.75 0]; 
legend({'Extracellular Autoxidation';'Gas';'Cellular'},'fontsize',16,'fontname','Helvetica')
ylabel('Total NO Consumption (\muM)','fontweight','bold')
xlabel('Time (hr)','fontweight','bold')
xlim([0 t(end)])
set(gca,'linewidth', 4,'fontsize',16,'fontname','Helvetica','box','off')

%Intracellular Cumulative NO Flux Plot
figure
A=area(t,cellular,'FaceColor','flat');
A(1).FaceColor=[0 1 0];
A(2).FaceColor=[0 0 0];
A(3).FaceColor=[0.5 0.5 0.5];
A(4).FaceColor =[1 1 0];
A(5).FaceColor =[1 0 1];
A(6).FaceColor =[0 1 1];


legend({'HMP';'NORV';'NRFA';'Iron-Sulfur Clusters';'Superoxide';'Cellular Autoxidation'},'fontsize',16,'fontname','Helvetica')
ylabel('Intracellular NO Consumption (\muM)','fontweight','bold')
xlabel('Time (hr)','fontweight','bold')
xlim([0 t(end)])
ylim([0 250])
set(gca,'linewidth', 4,'fontsize',16,'fontname','Helvetica','box','off')
