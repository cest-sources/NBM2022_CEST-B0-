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
B1= [0.75 2 4];
B1= [2];
B0=[3 7 9.4];
% T1= [0.939 1.222 1.429 	]  ;  % [Zhu 2014] http://hdl.handle.net/21.11116/0000-0001-32FE-9
T2= [0.062 0.037 0.029	]	;	 % [Zhu 2014] http://hdl.handle.net/21.11116/0000-0001-32FE-9

N= [1 80]; % number of pulses
for kk=1:numel(N)
    
for jj=1:numel(B1)
for ii=1:numel(B0)
seqid = sprintf('B0_%.1f_B1_%.2f_pulsed%d',B0(ii),B1(jj),N(kk));
seq_defs.B1rms=B1(jj);
seq_defs.n_pulses      = N(kk)            ; % number of pulses % 120 are better, 1 just for speedup test
seq_defs.shape ='sinc';
seq_defs.tp            = 15e-3            ; % pulse duration [s]
seq_defs.td            = 10e-3           ; % interpulse delay [s]
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
fa_sat        = seq_defs.B1rms*gyroRatio_rad*seq_defs.tp; % flip angle of sat pulse
if seq_defs.shape=="sinc"
    satPulse      = mr.makeSincPulse(fa_sat, 'Duration', seq_defs.tp, 'system', lims,'timeBwProduct', 2,'apodization', 0.15); % philips-like sinc
    [B1rms,B1cwae,B1cwae_pure,alpha]= calc_power_equivalents(satPulse,seq_defs.tp,seq_defs.td,0,gyroRatio_hz);
    satPulse.signal=satPulse.signal/B1rms*seq_defs.B1rms; % this line adjusts the pulse again to use power equivalent or B1rms metric
    [B1rms,B1cwae,B1cwae_pure,alpha]= calc_power_equivalents(satPulse,seq_defs.tp,seq_defs.td,0,gyroRatio_hz)
else
    satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', seq_defs.tp, 'system', lims);
end

satPulse      = resamplePulseForRLE(satPulse, 200); % resample pulse for reduced file size and io time

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
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated phase from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        
        % calc phase for next rf pulse
      % accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
       % accumPhase = mod(accumPhase + sign(currentOffset)*50*pi/180,2*pi); % 50° phase increment with offset polarity
       % accumPhase = mod(accumPhase + sign(currentOffset)*117*(np+1)*pi/180,2*pi); % quadr. 117° phase increment with offset polarity
        accumPhase = mod(accumPhase + sign(currentOffset)*113*pi/180,2*pi); % 113° phase increment with offset polarity:      works best!
      
       
        if np < seq_defs.n_pulses % delay between pulses
            if 1 % no gradspoil after each pulse
               seq.addBlock(mr.makeDelay(seq_defs.td)); % add delay
            else % gradspoil after each pulse
               seq.addBlock(mr.makeDelay(seq_defs.td-spoilDuration)); % add delay
               seq.addBlock(gxSpoil,gySpoil,gzSpoil);
            end
        end
    end
    % add the spoiling gradients
    seq.addBlock(gxSpoil,gySpoil,gzSpoil);
    % add the readout trigger event
    seq.addBlock(pseudoADC); 
end
toc
%% 
%%
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

%Psim=rmfield(Psim,'isochromats'); 
M_z{ii,jj} = simulate_pulseqcest(seq,Psim);
% generate Zref spectra of only water and MTm thus remove CEST pools from Psim
np = size(Psim.CESTPool,2)+1; Psim=rmfield(Psim,'CESTPool'); Psim.M=Psim.M([1 np+1 2*np+1, end],:); % manipulate the simulation parameters: remove all CEST pools
M_z_ref{ii,jj} = simulate_pulseqcest(seq,Psim);
end
end



%% plotting Z  and MTRLD
linest={'-','-','-'};
figure('Name',sprintf('Figure6_Part_%i_N=%i',kk,N(kk))),
subplot(4,1,1), 
for jj=numel(B1)
    set(gca,'ColorOrderIndex',1)
for ii=1:numel(B0)
w=seq_defs.offsets_ppm(2:end); % remove normalization scan
Z=M_z{ii,jj}(2:end)/M_z{ii,jj}(1);   % normalize by first scan
Z_ref=M_z_ref{ii,jj}(2:end)/M_z_ref{ii,jj}(1);   % normalize by first scan
% finally, plot the Z-spectrum 

hold on; grid on;
plot(w,Z,'LineStyle',linest{jj},'Displayname',sprintf('B_0=%.1f T , B_1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
% plot(w,Z_ref,'Displayname',sprintf('B0=%.1f T , B1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
xlabel('\Delta\omega [ppm]'); 
ylabel('Z(\Delta\omega)');
end

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
% finally, plot the MTR-spectrum 
subplot(4,1,jj+1), hold on; grid on;
plot(w,MTR,'LineStyle',linest{jj},'Displayname',sprintf('%.1f T',B0(ii))); set(gca,'xdir','reverse');
% plot(w,Z_ref,'Displayname',sprintf('B0=%.1f T , B1=%.2f µT',B0(ii),B1(jj))); set(gca,'xdir','reverse');
xlabel('\Delta\omega [ppm]');
ylabel('MTR_{LD}(\Delta\omega)');
ylim([0 0.09]);
end
l=legend('show'); l.FontSize = 8; l.Location='northeastoutside';
end

saveas(gcf, sprintf('Figure6_Part_%i_N=%i',kk,N(kk)));
end


%%
% uiopen('D:\root\LABLOG\FAU\26_NBM_CESTmultiB0\Fig_101\FIG_MTRLD_B0_v2.fig',1)
% uiopen('D:\root\LABLOG\FAU\26_NBM_CESTmultiB0\Fig_S1\data\FigS1_meas_B0_approx06µT.fig',1)


