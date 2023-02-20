function SelectMovingWindowSize(parameter, fileName)

SAD = zeros(10,2); %Sum of Absolute Differences Matrix

SAD(:,1) = [5;25;55;85;125;175;225;295;355;405];

for i = 1:10
    tmp_smoothed_data = sgolayfilt(parameter,2,SAD(i,1));
    SAD(i,2) = sum(abs(tmp_smoothed_data - parameter));
end

figure;

plot(SAD(:,1), SAD(:,2), '-ro','LineWidth',2);

xlabel('Window Size');
ylabel('SAD');

saveas(gcf, fileName, 'png');






    


