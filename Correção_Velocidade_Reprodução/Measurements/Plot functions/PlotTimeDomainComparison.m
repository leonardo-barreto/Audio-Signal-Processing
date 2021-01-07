function PlotTimeDomainComparison(time,refSig,newSigs)

for index = 1:length(newSigs)
    figure;
    plot(time,refSig);
    hold on;
    plot(time(1:length(newSigs{index})),newSigs{index});
    set(gca,'FontSize', 40)
    legend('Sinal original','Sinal revertido');
end
