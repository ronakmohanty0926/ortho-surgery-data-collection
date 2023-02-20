function plotData(testData, timeHaptics, timeForces, trueClass, predClass)

f1 = figure;
subplot(2,1,1)
hold on
for i = 1:length(testData)
    if trueClass(i,1) == 0
        scatter(timeHaptics(i,1), testData(i,1),'k');
    elseif trueClass(i,1) == 1
        scatter(timeHaptics(i,1), testData(i,1),'r');
    elseif trueClass(i,1) == 2
        scatter(timeHaptics(i,1), testData(i,1),'g');
    elseif trueClass(i,1) == 3
        scatter(timeHaptics(i,1), testData(i,1),'b');
    end
end

subplot(2,1,2)
for j = 1:length(testData)
    if predClass(j,1) == 0
        scatter(timeHaptics(j,1), testData(j,1),'k');
    elseif trueClass(j,1) == 1
        scatter(timeHaptics(i,1), testData(j,1),'r');
    elseif trueClass(j,1) == 2
        scatter(timeHaptics(i,1), testData(j,1),'g');
    elseif trueClass(j,1) == 3
        scatter(timeHaptics(j,1), testData(j,1),'b');
    end
end
hold off


