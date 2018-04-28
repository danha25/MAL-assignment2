Fs=8000; %Sampling frequency
nBits=16;
nChannels=1;
time=3; %Recording time in seconds

filePath='recordings/%s%d.mat';
word='dog'; % word you are recording
startIndex=1;   % start index of recording number 
repeat=3; % repeat recording


recObj = audiorecorder(Fs,nBits,nChannels);

recordingIndex=startIndex;
while recordingIndex<(startIndex+repeat)
    
    % Record audio signal
    disp('Start speaking.');
    recordblocking(recObj, time);
    disp('End of Recording.');

    % Play back the recording.
    play(recObj); % Listen if the recording is ok

    % Save the recording
    saveRecord = input('Do you want to save recording? Y/N [Y]: ', 's');
    if (strcmpi(saveRecord,'y') || isempty(saveRecord))
        % Store data in double-precision array.
        myRecording = getaudiodata(recObj);
        
        %Save recording to file. Filename format: "word""recording number".mat, ex. hello1.mat
        fileName=sprintf(filePath,word,recordingIndex);
        save(fileName,'myRecording') % Remember that the variable name will still be myRecording when you load the file
        disp(strcat(fileName, ' saved sucesfully'));
        recordingIndex=recordingIndex+1;
    end
end    
    
% Plot the recording in the time and frequency domains.
% L=Fs*time;
% fftMyRecording = fft(myRecording);
% psdTwoSided = abs(fftMyRecording/L);
% psdOneSided = psdTwoSided(1:L/2+1);
% psdOneSided(2:end-1) = psdOneSided(2:end-1)*2;
% 
% freqLabel = Fs*(0:(L/2))/L; %Frequency labels on the x-axis
% 
% figure()
% subplot(2,1,1)
% plot(myRecording);
% title('Time domain')
% subplot(2,1,2)
% plot(freqLabel,psdOneSided)
% title('Frequency domain')
