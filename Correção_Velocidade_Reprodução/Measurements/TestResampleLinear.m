function output = TestResampleLinear(chirp_original,Fs,windowSize)

numberWindows = ceil(length(chirp_original)/windowSize);
resample_factors = 1./linspace(1,2,numberWindows);



output = TimeVarying_ResampleV3(chirp_original,Fs,resample_factors,512,windowSize,0);

fileName = strcat('gymno_corrected_',num2str(windowSize),'.wav');

audiowrite(fileName,output,Fs);
end
