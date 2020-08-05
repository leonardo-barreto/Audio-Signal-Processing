function WriteAudioFiles(signalInformationStruct,inputName,fs)

    curveNames = signalInformationStruct.curveNames;

    for curveIndex = 1:length(signalInformationStruct.revertedSignals)

        currentsignal = signalInformationStruct.revertedSignals{curveIndex}./max(signalInformationStruct.revertedSignals{curveIndex});
        audioName = [inputName '_' curveNames{curveIndex} '.wav'];

        audiowrite(audioName,currentsignal,fs);

    end
    
end