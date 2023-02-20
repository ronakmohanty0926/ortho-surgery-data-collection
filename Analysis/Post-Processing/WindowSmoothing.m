function [signature,remPts,retPts]= WindowSmoothing(curve,winSize,nIterations)

signature = zeros(length(curve(:,1)),1);
smCurve = curve;
temp = curve;
for j = 1:1:nIterations
    for i = 2:1:length(curve(:,1))-1
        temp(i,:) = 0.5.*(smCurve(i-1,:)+smCurve(i+1,:));
    end
    
    for i = 1:1:length(curve(:,1))
        if(norm(smCurve(i,:)-temp(i,:)) > 0.001)
            signature(i) = signature(i)+1;
        end
    end
    smCurve = temp;
end

signature = signature./(max(signature)-min(signature));

remPts = [];
retPts = [];
for i = 1:1:length(curve(:,1))
    if(signature(i) < 0.2)
        retPts = [retPts;curve(i,:)];
    else
        remPts = [remPts;curve(i,:)];
    end
end

view(3);
plot3(curve(:,1), curve(:,2), curve(:,3),'r');
hold on;
plot3(smCurve(:,1), smCurve(:,2), smCurve(:,3),'b');
plot3(remPts(:,1), remPts(:,2), remPts(:,3),'ok');
plot3(retPts(:,1), retPts(:,2), retPts(:,3),'og');
daspect([1,1,1]);

figure;
plot(signature);

figure;
hist(signature);