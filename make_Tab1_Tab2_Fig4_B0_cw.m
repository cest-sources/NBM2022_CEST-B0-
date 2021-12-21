%% Pulseq-CEST
% The following live script demonstrates how to generate an pulseq-file for 
% an APT-weighted CEST experiment with Pulseq-CEST.
%% Installation
% The Pulseq-CEST MATLAB Code can be cloned from git and is set up using the 
% Installation script

% if you do not have git installed you can use the commented lines to unzip
% the code directly from GitHub:
% unzip("https://github.com/kherz/pulseq-cest/archive/master.zip");
% movefile('pulseq-cest-master', 'pulseq-cest');
 if exist('pulseq-cest', 'dir')
        disp('pulseq-cest-library already installed, skip...')
 else
    system('git clone -b v1.0.0 https://github.com/kherz/pulseq-cest'); 
    cd pulseq-cest;
    install_pulseqcest;
    cd ..
 end
 
pulseqCEST_simlib=[pwd '\pulseq-cest-library\sim-library\'];


%% 
% The _seq_def_ struct contains everything that gets later stored as a definition 
% in the pulseq-file. This can contain anything you would need later for post-processing 
% etc.
%%
clear M_z
B1= [0.75 2 4];
B0=[3 7 9.4];
T1_zhu= [0.939 1.222 1.429 	]      ;                 % [Zhu 2014] http://hdl.handle.net/21.11116/0000-0001-32FE-9
T2= [0.062 0.037 0.029	]	;	% [Zhu 2014] http://hdl.handle.net/21.11116/0000-0001-32FE-9

R1A=0.4; R1C=12.2./B0;
fc=0.289; kac = 1.38; kca= 1.38/fc;
R1obs_wang_approx = (R1A + fc*R1C )/(1 + fc);
R1obs_wang =    0.5*( kac + kca + R1A+ R1C - sqrt(( kac + kca + R1A + R1C ).^2 - 4*( kca*R1A + kac.*R1C + R1A.*R1C )));
1./R1obs_wang
T1_zhu
for jj=1:numel(B1)

for ii=1:3
seqid = sprintf('B0_%.1f_B1_%.2f_cw',B0(ii),B1(jj));

seq_defs.B1=B1(jj);
seq_defs.n_pulses      = 1             ; % number of pulses
seq_defs.shape ='block'
seq_defs.tp            = 2000e-3           ; % pulse duration [s]
seq_defs.td            = 0e-3            ; % interpulse delay [s]
seq_defs.Trec          = 3.5             ; % recovery time [s]
seq_defs.Trec_M0       = 3.5             ; % recovery time before M0 [s]
seq_defs.M0_offset     = -1560           ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -100 -50 -30 -20 -10 -8:0.1:8 10 20 30 50 100]; % offset vector [ppm]
seq_defs.offsets_ppm   = [seq_defs.M0_offset -8:0.1:8 ]; % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
                         seq_defs.td ;  % saturation time [s]
seq_defs.B0            = B0(ii)               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id

% init sequence
% seq = mr.Sequence();
seq = SequenceSBB();
%% set definitions in the seq file
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end

lims = getScannerLimits();
%% 
% In the next step we create all blocks needed: the saturation pulses, spoiler 
% gradients, and ADC events. Pulseq comes with some basic pulse shapes that can 
% be further adapted with additional parameters. For this example we use a long 
% block pulse which is also called contiouus wave irradiaton. The spoiler gradients 
% are strong trapezoidial gradients in all directions. The pseudo-ADC event defines 
% the time pint at which the signal is sampled in the simulation, or the readout 
% sequence is called on the scanner.

% rf pulses
gyroRatio_hz  = 42.5764;               % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;     % [rad/uT]
fa_sat        = seq_defs.B1*gyroRatio_rad*seq_defs.tp; % flip angle of sat pulse
if seq_defs.shape=="sinc"
    satPulse      = mr.makeSincPulse(fa_sat, 'Duration', seq_defs.tp, 'system', lims,'timeBwProduct', 4,'apodization', 0.15); % philips-like sinc
else
    satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', seq_defs.tp, 'system', lims);
end
% spoilers
spoilRiseTime = 1e-3;
spoilDuration = 4500e-6 + spoilRiseTime; % [s]
[gxSpoil, gySpoil, gzSpoil] = makeSpoilerGradients(lims, spoilDuration, spoilRiseTime);

% pseudo adc, not played out
pseudoADC = mr.makeAdc(1,'Duration', 1e-3);
%% 
% Now we have all the objects we need and can fill the pulseq-file by looping 
% through the offsets we want to simulate / measure.
%%
tic

% pulseq uses offsets in Hz
offsets_Hz = seq_defs.offsets_ppm*gyroRatio_hz*seq_defs.B0;

% loop through offsets
for currentOffset = offsets_Hz
    % set frequency offset of the pulse
    satPulse.freqOffset = currentOffset;  
    accumPhase=0;
    for np=1:seq_defs.n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        
        seq.addBlock(satPulse) % add sat pulse
        
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < seq_defs.n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(seq_defs.td)); % add delay
        end
    end
    % add the spoiling gradients
    seq.addBlock(gxSpoil,gySpoil,gzSpoil);
    % add the readout trigger event
    seq.addBlock(pseudoADC); 
end
toc
%% 
% Psim = readSimulationParameters('./WM_3T_001_bmsim.yaml');
% Psim = readSimulationParameters('./pulseq-cest-library\sim-library\WM_3T_Liu13.yaml');

if 0 % orig
Psim = readSimulationParameters([pulseqCEST_simlib 'WM_3T_Stanisz2005_5pool_bmsim.yaml']);
Psim.WaterPool.R1=1/T1(ii);
Psim.WaterPool.R2=1/T2(ii);
else % wang
Psim = readSimulationParameters([pulseqCEST_simlib 'WM_3T_Wang2020_5pool_bmsim.yaml']);
Psim.WaterPool.R1=0.4;
Psim.WaterPool.R2=1/T2(ii);
Psim.MTPool.R1=12.2/B0(ii);
% Psim.MTPool.R2=Psim.MTPool.R2 * T2(ii)/T2(1); warning('MT T2 tweak');
end

Psim.Scanner.B0=seq_defs.B0; % set same B0 as used in sequence definition

% Psim.MTPool.R1=12.2/B0(ii);
% Psim.MTPool.R2=3/B0(ii);
% Psim.MTPool.R2=Psim.MTPool.R2 * T2(1)/T2(ii); warning('MT T2 tweak');
% F(1)=Psim.MTPool.f ; F(2)=Psim.MTPool.f * 1.5200; F(3)=Psim.MTPool.f * 1.5200*9.4/7; 
% Psim.MTPool.f= F(ii);

M_z{ii,jj} = simulate_pulseqcest(seq,Psim);
np = size(Psim.CESTPool,2)+1; Psim=rmfield(Psim,'CESTPool'); Psim.M=Psim.M([1 np+1 2*np+1, end],:); % manipulate the simulation parameters: remove all CEST pools
M_z_ref{ii,jj} = simulate_pulseqcest(seq,Psim);
%% 
% The simulation returns the Z-magnetization which we can now plot together 
% with the MTRasym.
end 
end

%% plotting Z  and MTRLD
linest={'-','-','-'};
figure('Name','Figure 4: Z-spectra and MTRLD simulation -B0,B1'),
subplot(4,1,1), 
for jj=1:numel(B1)
    set(gca,'ColorOrderIndex',1)
for ii=1:numel(B0)
w=seq_defs.offsets_ppm(2:end); % remove normalization scan
Z=M_z{ii,jj}(2:end)/M_z{ii,jj}(1);   % normalize by first scan
Z_ref=M_z_ref{ii,jj}(2:end)/M_z_ref{ii,jj}(1);   % normalize by first scan
% finally, plot the Z-spectrum 

hold on; grid on;
set(gca,'ColorOrderIndex',ii)
plot(w,Z,'LineStyle',linest{jj},'Displayname',sprintf('B_0=%.1f T , B_1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
% plot(w,Z_ref,'Displayname',sprintf('B0=%.1f T , B1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
if 0 % plot MTRasym
set(gca,'ColorOrderIndex',ii)
plot(w,Z(end:-1:1)-Z,'LineStyle',linest{jj},'Displayname',sprintf('B_0=%.1f T , B_1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
end
xlabel('\Delta\omega [ppm]'); 
ylabel('Z(\Delta\omega)');
end
xlim([-8  8]);

end
legend({'3.0 T','7.0 T','9.4 T'},'FontSize',8, 'Location','northeastoutside');
set(gcf,'Position',[587   262   560   720]);
t = annotation('textbox','String','a)','Position',[0 0.85 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','b)','Position',[0 0.64 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','c)','Position',[0 0.42 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
t = annotation('textbox','String','d)','Position',[0 0.2 0.1 0.1]);
t.FontSize = 14; t.LineStyle='None';
set(groot,'defaultLineLineWidth',1.0)

% plotting MTRLD
for jj=1:numel(B1)
    t = annotation('textbox','String',sprintf('%.2f µT',B1(jj)),'Position',[0.63 0.83-0.22*jj 0.1 0.1]); t.LineStyle='None';
    t = annotation('textbox','String',sprintf('%.2f µT',B1(jj)),'Position',[0.64 0.845-0.03*jj 0.1 0.1]); t.LineStyle='None'; t.FontSize = 8;
for ii=1:numel(B0)
w=seq_defs.offsets_ppm(2:end); % remove normalization scan
Z=M_z{ii,jj}(2:end)/M_z{ii,jj}(1);   % normalize by first scan
Z_ref=M_z_ref{ii,jj}(2:end)/M_z_ref{ii,jj}(1);   % normalize by first scan
MTR=Z_ref-Z;
% finally, plot the Z-spectrum 
subplot(4,1,jj+1), hold on; grid on;
plot(w,MTR,'LineStyle',linest{jj},'Displayname',sprintf('%.1f T',B0(ii))); set(gca,'xdir','reverse');
% plot(w,Z_ref,'Displayname',sprintf('B0=%.1f T , B1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
xlabel('\Delta\omega [ppm]');
ylabel('MTR_{LD}(\Delta\omega)');
ylim([0 0.09]); xlim([-8 8]);
end
l=legend('show'); l.FontSize = 8; l.Location='northeastoutside';
end

saveas(gcf,'Figure4_Z-MTRLD_B0_B1');

%% comparison to real data - only works if you are MZ
% uiopen('D:\root\LABLOG\FAU\26_NBM_CESTmultiB0\Fig_101\FIG_MTRLD_B0_v2.fig',1)
% uiopen('D:\root\LABLOG\FAU\26_NBM_CESTmultiB0\Fig_S1\data\FigS1_meas_B0_approx06µT.fig',1)

%% Table 1 - overview of simulation parameters
Psim = readSimulationParameters([pulseqCEST_simlib 'WM_3T_Wang2020_5pool_bmsim.yaml']);
clear T
T=table();
W = Psim.WaterPool;
MT= Psim.MTPool;
T(1,:) = {W.R1 W.R2 W.f 0 0}
T(2,:) ={MT.R1 MT.R2 MT.f MT.dw MT.k}
T.Properties.RowNames={'Water' 'ssMT'};
for ii=1:numel(Psim.CESTPool)
T(ii+2,:) ={Psim.CESTPool(ii).R1 Psim.CESTPool(ii).R2 Psim.CESTPool(ii).f Psim.CESTPool(ii).dw Psim.CESTPool(ii).k}
T.Properties.RowNames{ii+2}=Psim.CESTPool(ii).id;
end

if 1 % use relax times instead of rates
T.Properties.VariableNames{1}='T1'; T.Properties.VariableNames{2}='T2';
T{:,1}=1./T{:,1}; T{:,2}=1./T{:,2}; 
end
T = T(:,[4 5 3 1 2])
writetable(T,'Table1_Wang.xls','WriteRowNames',true);


%% Table 2
B1_idx=2; CEST_idx=2;   % 2µT , guanidine

R1A=0.4; R1C=12.2./B0;
fc=0.289; kac = 1.38; kca= 1.38/fc;
R1obs_wang_approx = (R1A + fc*R1C )/(1 + fc);
R1obs_wang =    0.5*( kac + kca + R1A+ R1C - sqrt(( kac + kca + R1A + R1C ).^2 - 4*( kca*R1A + kac.*R1C + R1A.*R1C )));
T1obs=1./R1obs_wang;

Psim = readSimulationParameters([pulseqCEST_simlib 'WM_3T_Wang2020_5pool_bmsim.yaml']);
C=Psim.CESTPool(CEST_idx);  % choose the CEST pool : 1 amide, 2 guanidine 
w1=gyroRatio_rad*B1(B1_idx);
alpha = w1^2./(C.k*(C.k+C.R2)+w1^2)
C.k*alpha
clear T
T=table();
T{:,1} = {sprintf('%s: %.1fppm k=%.1fHz B1=%.2fµT alpha=%.1f%%',C.id, C.dw,C.k,B1(B1_idx),alpha*100),'',''}'; T.Properties.VariableNames{1}='Table2';
T{:,2} = {'','',''}'; T.Properties.VariableNames{2}='at';
T{:,3} = B0'; T.Properties.VariableNames{3}='B0';
T{:,4} = 1./R1A; T.Properties.VariableNames{4}='T1w';
T{:,5} = 1./R1C'; T.Properties.VariableNames{5}='T1mt';
T{:,6} = T1obs'; T.Properties.VariableNames{6}='T1obs';
T{:,7} = T2'; T.Properties.VariableNames{7}='T2';
[~, idx]=find(seq_defs.offsets_ppm>=C.dw); idx=idx(1); seq_defs.offsets_ppm(idx)
T{:,8} = [ M_z_ref{1,2}(idx) M_z_ref{2,2}(idx) M_z_ref{3,2}(idx)]'.^2; T.Properties.VariableNames{8}='sigma';
T{:,9}=T{:,6}.*T{:,8}.*alpha.*C.k*C.f*100 ;   T.Properties.VariableNames{9}='MTRLD';
T{:,10}=T{:,6}.*T{:,8}.*alpha.*C.k ;    T.Properties.VariableNames{10}='PTE';
T
writetable(T,'Table2_Wang.xls','WriteRowNames',true);

