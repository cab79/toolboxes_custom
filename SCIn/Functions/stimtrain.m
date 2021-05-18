function h = stimtrain(h,opt)

if ~isfield(h,'stim')
    h.stim = struct;
end

switch h.Settings.stim(h.sn).control

    case {'labjack','LJTick-DAQ','T7'}
        
        switch opt
            case 'create'
                h=trial(h,'set&construct');
                return
        end
        
%         if ~isfield(h,'ljHandle')
%             try
%                 h.ljHandle = get(h.ljhandle, 'Value');
%             end
%         end
        
        if h.Settings.labjack
            if ~isunix
                if strcmp(h.Settings.stim(h.sn).control,{'labjack','LJTick-DAQ'})
                    ljud_LoadDriver; % Loads LabJack UD Function Library
                    ljud_Constants; % Loads LabJack UD constant file
                elseif strcmp(h.Settings.stim(h.sn).control,{'T7'})
                    t = h.ljmAsm.AssemblyHandle.GetType('LabJack.LJM+CONSTANTS'); % Creating an object to nested class LabJack.LJM.CONSTANTS
                    LJM_CONSTANTS = System.Activator.CreateInstance(t);
                end
            elseif isempty(h.ljHandle.cal)
                getCal(h.ljHandle);
            end
        end
        
        if h.Settings.stim(h.sn).chanforLJ
            port = h.Settings.stim(h.sn).chan;
        else
            port=6;
        end
        
        switch opt
            case 'setup'
                % if GUI is set to intensity > 0, check this is not
                % accidental
                if isfield(h,'inten_mean_gui')
                    inten_mean = str2double(h.inten_mean_gui);
                    if inten_mean>0
                        choice = questdlg(['Is ' num2str(inten_mean) ' the correct starting intensity?'], ...
                        'Intensity', ...
                        'Yes','No','Yes');
                        % Handle response
                        switch choice
                            case 'No'
                                dbclear all
                                error('Quitting: wrong intensity set in GUI')
                        end
                    end
                end
            
            case 'create'
               h=trial(h,'set&construct');
                
            case 'calc'
                h=trial(h,'set');
                
            case 'getsample'
                h.currentsample = round((h.ct-h.out.stimtime{1})*h.Settings.fs);
            
            case 'set'
                
                % if selecting the intensity from a time sample of a
                % waveform
                if strcmp(h.Settings.design,'continuous')
                    if isfield(h,'currentsample')
                        if h.currentsample>0
                            h.stim(h.sn).inten = h.Seq.stimseq(h.currentsample);
                            h.out.stimseq_record(h.currentsample) = h.stim(h.sn).inten;
                        else
                            h.stim(h.sn).inten = 0;
                        end
                    else
                        h.stim(h.sn).inten = 0;
                    end
                end
                
                
                %Set DACA 
                if h.Settings.labjack
                    if isunix
                        h.ljHandle.setDAC(h.stim(h.sn).inten*h.Settings.DAC_multiply);
                        WaitSecs(0.1);
                    else
                        if strcmp(h.Settings.stim(h.sn).control,'T7')
                            LabJack.LJM.eWriteName(h.ljHandle, 'DAC1', h.stim(h.sn).inten*h.Settings.DAC_multiply);
                        elseif strcmp(h.Settings.stim(h.sn).control,'LJTick-DAQ')
                            try
                                Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, h.stim(h.sn).inten*h.Settings.DAC_multiply, 0); 
                                if Error~=0
                                    error('error')
                                end
                            catch% to use DAC0 port
                                Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, h.stim(h.sn).inten*h.Settings.DAC_multiply, 0,0);
                                Error_Message(Error)
                                Error = ljud_GoOne(h.ljHandle);
                                Error_Message(Error)
                            end
                        elseif strcmp(h.Settings.stim(h.sn).control,'labjack')
                            Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, h.stim(h.sn).inten*h.Settings.DAC_multiply, 0,0);
                            Error_Message(Error)
                            Error = ljud_GoOne(h.ljHandle);
                            Error_Message(Error)
                        end
                    end
                end

                %voltage = 80;
                %current = h.stim(h.sn).inten*100;
                %time = 50/1e6 * h.Settings.npulses_train;
                %mJ = voltage * current * time
               
                
            case 'start'
                
                if ~h.Settings.labjack
                    return
                end

                if h.Settings.stim(h.sn).labjack_timer && ~isempty(h.Settings.stim(h.sn).p_freq)
                    if h.Settings.stim(h.sn).p_freq>=h.LJfreqtable(1,1)
                        %Set the timer/counter pin offset to 6, which will put the first timer/counter on FIO6. 
                        Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                        Error_Message(Error)

                        % get row of table. Columns are Hz, base clock, clock divisor, and timer value. 
                        r = dsearchn(h.LJfreqtable(:,1),h.Settings.stim(h.sn).p_freq);
                        % get parameters
                        freq = h.LJfreqtable(r,1);
                        base = h.LJfreqtable(r,2);
                        div = h.LJfreqtable(r,3);
                        val = h.LJfreqtable(r,4);
                        disp(['actual freq is ' num2str(freq)]);
                        
                        if base == 1e6; basenum = 23;%   //1 MHz clock base w/ divisor (no Counter0)
                        elseif base == 4e6; basenum = 24;%   //4 MHz clock base w/ divisor (no Counter0)
                        elseif base == 12e6; basenum = 25;%  //12 MHz clock base w/ divisor (no Counter0)
                        elseif base == 48e6; basenum = 26;%  //48 MHz clock base w/ divisor (no Counter0)
                        end

                        %use 48MHz clock base with divisor = 48 to get 1 MHz timer clock: 
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_BASE, basenum, 0, 0);
                        Error_Message(Error)
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_DIVISOR, div, 0, 0); 
                        Error_Message(Error)

                        %Enable 2 timers.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chNUMBER_TIMERS_ENABLED, 2, 0, 0); 
                        Error_Message(Error)

                        %Configure Timer0 as Frequency out.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 0, LJ_tmFREQOUT, 0, 0); 
                        Error_Message(Error)

                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 0, val, 0, 0); 
                        Error_Message(Error)
                        
                    else % use less accurate method for low freq
                        %Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, 0.5, 0,0);
                        %Error_Message(Error)
                        %Error = ljud_GoOne(h.ljHandle);
                        %Error_Message(Error)
                    
                        %Set the timer/counter pin offset to 6, which will put the first timer/counter on FIO6. 
                        Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                        Error_Message(Error)

                        %Enable 2 timers.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chNUMBER_TIMERS_ENABLED, 2, 0, 0); 
                        Error_Message(Error)

                        % use 12MHz clock base
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_BASE, LJ_tc12MHZ_DIV, 0, 0);
                        Error_Message(Error)

                        % control freq output with a divisor from 6 to 183 to get an output frequency of about 30.5Hz to 1Hz respectively: 
                        div = round(12e6 / (h.Settings.stim(h.sn).p_freq * 2^16));
                        freq = 12e6/div/2^16;
                        disp(['actual freq is ' num2str(freq)]);
                        %div = round(12e6 / (0.1 * 2^16));
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_CONFIG, LJ_chTIMER_CLOCK_DIVISOR, div, 0, 0); 
                        Error_Message(Error)

                        %Configure Timer0 as 16-bit PWM.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 0, LJ_tmPWM16, 0, 0); 
                        Error_Message(Error)

                        %Initialize the 16-bit PWM with a 50% duty cycle.
                        Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 0, 32767, 0, 0); 
                        Error_Message(Error)
                    end

                    %Configure Timer1 as timer stop:
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_MODE, 1, LJ_tmTIMERSTOP, 0, 0);
                    Error_Message(Error)

                    %set number of pulses: 
                    Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_TIMER_VALUE, 1, h.Settings.stim(h.sn).npulses_train, 0, 0);
                    Error_Message(Error)

                    %Execute the requests. 
                    Error = ljud_GoOne(h.ljHandle);
                    Error_Message(Error)
                    disp('running')
                    
                elseif isfield(h.Settings.stim(h.sn),'wavetype') && strcmp(h.Settings.stim(h.sn).wavetype,'digital')
                    % assumes T7 control using StreamOut
                    
%                     % test
%                     value = 5;  % 2.5 V
%                     LabJack.LJM.eWriteName(h.ljHandle, name, value);
%                     pause(0.1)

                    % seems to need zeroing
                    LabJack.LJM.eWriteName(h.ljHandle, port, 0);
                    
                    % Setup stream-out
                    numAddressesOut = 1;
                    aNamesOut = NET.createArray('System.String', numAddressesOut);
                    aNamesOut(1) = port;
                    aAddressesOut = NET.createArray('System.Int32', numAddressesOut);
                    aTypesOut = NET.createArray('System.Int32', numAddressesOut);  % Dummy
                    LabJack.LJM.NamesToAddresses(numAddressesOut, aNamesOut, ...
                        aAddressesOut, aTypesOut);

                    % Allocate memory for the stream-out buffer
                    max_buffer = 16384;
                    if length(h.mwav)>(max_buffer/2-1)
                        warning('stimtrain.m: T7 buffer exceeds maximum - truncating')
                        h.mwav = h.mwav(1:max_buffer/2-1);
                    end
                    buffer_size = max(32,pow2(ceil(log2(length(h.mwav))))*2); % min of 32
                    LabJack.LJM.eWriteName(h.ljHandle, 'STREAM_OUT0_ENABLE', 0);
                    LabJack.LJM.eWriteName(h.ljHandle, 'STREAM_OUT0_TARGET', aAddressesOut(1));
                    LabJack.LJM.eWriteName(h.ljHandle, 'STREAM_OUT0_BUFFER_SIZE', buffer_size);
                    LabJack.LJM.eWriteName(h.ljHandle, 'STREAM_OUT0_ENABLE', 1);
                    LabJack.LJM.eWriteName(h.ljHandle, 'STREAM_OUT0_LOOP_NUM_VALUES', 0);
                    LabJack.LJM.eWriteName(h.ljHandle, 'STREAM_OUT0_SET_LOOP', 1);

                    % Load the waveform data points
                    LabJack.LJM.eWriteNameArray(h.ljHandle, 'STREAM_OUT0_BUFFER_F32', min(length(h.mwav),max_buffer/2), h.mwav(1:min(length(h.mwav),max_buffer)), 0);
                    
                    % Setup stream-out
                    numAddressesOut = 1;
                    aNamesOut = NET.createArray('System.String', numAddressesOut);
                    aNamesOut(1) = 'STREAM_OUT0';
                    aAddressesOut = NET.createArray('System.Int32', numAddressesOut);
                    aTypesOut = NET.createArray('System.Int32', numAddressesOut);  % Dummy
                        LabJack.LJM.NamesToAddresses(numAddressesOut, aNamesOut, ...
                            aAddressesOut, aTypesOut);

                    % Configure and start stream
                    if (h.freq>100000)
                        error('reduce scan rate to <100,000')
                    end
                    aScanList = aAddressesOut(1); % FIO0
                    scanRate = h.freq; % Hz
                    scansPerRead = scanRate/2; % eStreamRead frequency = ScanRate / ScansPerRead
                    numAddresses = length(aScanList);

                    [~, scanRate] = LabJack.LJM.eStreamStart(h.ljHandle, scansPerRead, ...
                        numAddresses, aScanList, scanRate);
                    disp(['Stream started with a scan rate of ' num2str(scanRate) ...
                                  ' Hz.'])

%                     % optional plotting
%                     if isfield(h,'fig'); close(h.fig); end
%                     incr=1/(length(h.mwav));
%                     x = (incr:incr:length(h.mwav)*incr)*(length(h.mwav)/scanRate);
%                     h.fig=figure; plot(x,h.mwav);
%                     xlabel('seconds')
%                     title('stimulus train waveform')
                    
                    WaitSecs(length(h.mwav)/scanRate)

                    % stop stream
                    disp('Stop Stream')
                    LabJack.LJM.eStreamStop(h.ljHandle);
                else
                    % pulse train instruction
                    if isempty(h.Settings.stim(h.sn).p_freq) || h.Settings.stim(h.sn).p_freq == 0
                       p_freq = 10000; % 0.1 ms delay (DS8R can detect down to 0.01ms)
                    else
                       p_freq = h.Settings.stim(h.sn).p_freq; 
                    end
                    if isunix
                        t1=GetSecs;
                        for pr = 1:h.Settings.stim(h.sn).npulses_train % train
                            
                            h.ljHandle.timedTTL(port,((1000000/p_freq)/2));
                            
                        end
                        t2=GetSecs;
                        disp(num2str(t2-t1))
                    else
                        for pr = 1:h.Settings.stim(h.sn).npulses_train % train
                            Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,port,1,0,0); % 
                            Error_Message(Error)

                            Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,port,round((1000000/p_freq)/2),0,0); % Actual resolution is 64 microseconds.
                            Error_Message(Error)

                            Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_DIGITAL_BIT,port,0,0,0);
                            Error_Message(Error)

                            Error = ljud_AddRequest(h.ljHandle,LJ_ioPUT_WAIT,port,round((1000000/p_freq)/2),0,0); % Actual resolution is 64 microseconds.
                            Error_Message(Error)
                        end
                        %Execute the stimulus train
                        t1=GetSecs;
                        Error = ljud_GoOne(h.ljHandle);
                        Error_Message(Error)
                        %ljud_GetResult(ljHandle, LJ_ioGET_DIGITAL_BIT, 7, @Value)
                        t2=GetSecs;
                    end
                    disp(['pulses per stim: ' num2str(pr) ';  labjack stim length: ' num2str(t2-t1)]);
                end
                
            case 'stop'
                
                if isunix
                    h.ljHandle.setFIO(0,0)
                else
                    if strcmp(h.Settings.stim(h.sn).control,'LJTick-DAQ')
                        Error = ljud_ePut(h.ljHandle, LJ_ioTDAC_COMMUNICATION, LJ_chTDAC_UPDATE_DACA, 0, 0); 
                        Error_Message(Error)
                    elseif strcmp(h.Settings.stim(h.sn).control,'labjack')
                        %try
                            Error = ljud_AddRequest(h.ljHandle, LJ_ioPUT_DAC, 0, 0, 0,0);
                            Error_Message(Error)
                            Error = ljud_GoOne(h.ljHandle);
                            Error_Message(Error)
                        %end
                        if isfield(h.Settings.stim(h.sn),'labjack_timer')
                            if h.Settings.stim(h.sn).labjack_timer
                                Error = ljud_AddRequest (h.ljHandle,  LJ_ioPUT_CONFIG, LJ_chTIMER_COUNTER_PIN_OFFSET, port, 0, 0);
                                Error_Message(Error)
                                %Execute the requests. 
                                Error = ljud_GoOne(h.ljHandle);
                                Error_Message(Error)
                            end
                        end
                    elseif strcmp(h.Settings.stim(h.sn).control,'T7')
                        Error = LabJack.LJM.eWriteName(h.ljHandle, 'DAC0', 0);
                        LabJack.LJM.Close(h.ljHandle);
                    end

                    
                end
                
        end
         
    case 'serial'
        switch opt
            case 'setup'
                opt = 'spt1';
                open_serial(h,opt);
                opt = 'spt';
                open_serial(h,opt);
            case 'set'
                h=trial(h,'set');
                global spt
                if (h.stim(h.sn).inten~=0 && h.stim(h.sn).inten<=255)
                    fprintf(spt,'%s', h.stim(h.sn).inten);
                else
                    error('invalid Intensity level')
                end

            case 'start'
                % trigger stimulator and mark EEG at the same time
                global spt1
                TriggerNum = h.Seq.signal(h.sn,h.i);
                fprintf(spt1,num2str(32+TriggerNum));
        end

    case 'audioplayer'
        if ~exist('opt','var')
            opt = 'run';
        end
        
        switch opt
            case 'run'
                h=trial(h,'set&construct');
                h.Seq.aud = audioplayer(h.Seq.stimseq', h.Settings.fs);
                play(h.Seq.aud);
                %sound(h.Seq.stimseq', h.Settings.fs, 16);
                %pause(0.2)
            
            case 'create'
                h=trial(h,'set&construct');
                h.Seq.aud = audioplayer(h.Seq.stimseq', h.Settings.fs);
            
            case 'start'
                play(h.Seq.aud);
                h.playstart = GetSecs;
                
            case 'pause'
                pause(h.Seq.aud);

            case 'resume'
                resume(h.Seq.aud);

            case 'stop'
                stop(h.Seq.aud);
                
            case 'getsample'
                h.currentsample=get(h.Seq.aud,'CurrentSample');
                h.totalsamples=get(h.Seq.aud,'TotalSamples');
        end
        
    case 'PsychPortAudio'
        switch opt
            case 'setup'
                h = PTBaudio(h);
                
            case 'getsample'
                s = PsychPortAudio('GetStatus', h.pahandle);
                if s.Active == 1
                    h.currentsample=s.ElapsedOutSamples;
                    %h.totalsamples=;
                end
                
            case 'create' 
                h=trial(h,'set&construct');
                %if strcmp(h.Settings.design,'trials')
                %    PsychPortAudio('FillBuffer', h.pahandle, h.Seq.stimseq);
                %    
                %elseif strcmp(h.Settings.design,'continuous')
                %    h.pabuffer = PsychPortAudio('CreateBuffer', h.pahandle, h.Seq.stimseq);% Engine still running on a schedule?
                %   
                %end
                
            case 'start' 
                if strcmp(h.Settings.design,'trials')
                    PsychPortAudio('FillBuffer', h.pahandle, h.Seq.stimseq);
                    PsychPortAudio('Start', h.pahandle, 1, 0, 1);
                    
%                     % optional plotting
%                     if isfield(h,'fig'); close(h.fig); end
%                     incr=1/(size(h.Seq.stimseq,2));
%                     x = (incr:incr:size(h.Seq.stimseq,2)*incr)*(size(h.Seq.stimseq,2)/h.Settings.fs);
%                     h.fig=figure; plot(x,h.Seq.stimseq(1,:));
%                     xlabel('seconds')
%                     title('stimulus train waveform')

                elseif strcmp(h.Settings.design,'continuous')
                    h.pabuffer = PsychPortAudio('CreateBuffer', h.pahandle, h.Seq.stimseq);
                   
                    s = PsychPortAudio('GetStatus', h.pahandle);
                    if s.Active == 0 %&& ~isfield(h,'i') % new run
                        PsychPortAudio('UseSchedule', h.pahandle, 1, size(h.Seq.signal,2));
                        PsychPortAudio('AddToSchedule', h.pahandle, h.pabuffer);
                        h.playstart = PsychPortAudio('Start', h.pahandle, 0, 0, 1);
                    %elseif s.Active == 0 
                    %    error('increase ntrialsahead in Settings')
                    else
                        PsychPortAudio('AddToSchedule', h.pahandle, h.pabuffer);
                        disp('new trial(s) added to schedule')
                    end
                end
        end
    case 'ptb_visual'
        switch opt
            case 'setup'
                h = PTBvisual(h);

            case 'getsample'
                % not needed

            case 'create' 
                h=trial(h,'set');

            case 'start' 
               h=VisualStim(h,'stim');
        end
end
end

