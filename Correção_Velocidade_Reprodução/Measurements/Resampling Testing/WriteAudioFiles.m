function WriteAudioFiles(signalInformationStruct,inputName,fs)

    curveNames = signalInformationStruct.curveNames;

    for curveIndex = 1:length(signalInformationStruct.revertedSignals)

        distortedSignal = signalInformationStruct.distortedSignals{curveIndex}./max(signalInformationStruct.distortedSignals{curveIndex});
        audioName = [inputName '_dist_' curveNames{curveIndex} '.wav'];

        audiowrite(audioName,distortedSignal,fs);

        revertedSignal = signalInformationStruct.revertedSignals{curveIndex}./max(signalInformationStruct.revertedSignals{curveIndex});
        audioName = [inputName '_rev_' curveNames{curveIndex} '.wav'];

        audiowrite(audioName,revertedSignal,fs);

    end
    
end