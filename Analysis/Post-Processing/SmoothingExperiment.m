% function [signature,remPts,retPts]= SmoothingExperiment(curve,nIterations)
% function SmoothingExperiment(curve,nIterations)
% function [arcLength, signature] = SmoothingExperiment(curve, drill_states)
function [arcLength, signature] = SmoothingExperiment(curve)

% dataLen = length(curve(:,1));
% dataLen = 1200;
% [curve, drill_states] = resample(ogcurve,ogdrill_states,dataLen);

dist = zeros(length(curve(:,2)),1);
for k = 2:length(curve(:,2))
    dist(k,1) = norm(curve(k,:) - curve(k-1,:));
end

arcLength = zeros(length(curve(:,2)),1);
S = sum(dist);
tot_dist = dist(1,1);
for m = 2:length(dist)
    tot_dist = tot_dist + dist(m,1);
    arcLength(m,1) = tot_dist./S;
end

signature = zeros(length(curve(:,1)),1);
smCurve = curve;
temp = curve;
while(1)    
    for i = 2:1:length(curve(:,1))-1
        temp(i,:) = 0.5.*(smCurve(i-1,:)+smCurve(i+1,:));
    end    
    globalCount = 0; 
    for i = 1:1:length(curve(:,1))
        if(norm(smCurve(i,:)-temp(i,:)) > 0.000001)
%         if(norm(smCurve(i,:)-temp(i,:)) > 0.000001)
            globalCount = globalCount + 1;
            signature(i) = signature(i)+1;
        end
    end
    smCurve = temp;
    
    if (globalCount == 0)
        break;
    end
end

% for m = 1:1:1000  
%     for i = 2:1:length(curve(:,1))-1
%         temp(i,:) = 0.5.*(smCurve(i-1,:)+smCurve(i+1,:));
%     end    
% %     globalCount = 0; 
%     for i = 1:1:length(curve(:,1))
%         if(norm(smCurve(i,:)-temp(i,:)) > 0.00001)
% %         if(norm(smCurve(i,:)-temp(i,:)) > 0.000001)
% %             globalCount = globalCount + 1;
%             signature(i) = signature(i)+1;
%         end
%     end
%     smCurve = temp;
%     
% %     if (globalCount == 0)
% %         break;
% %     end
% end


% Normalization
signature = (signature -  min(signature))./(max(signature)-min(signature));

% cort1count = 0;
% cort2count = 0;
% for i = 1:1:length(curve(:,1))
%     if(drill_states(i,1) == 1)
%         cort1count = cort1count + 1;
%         cort1sig(cort1count, 2) = signature(i,1);
%         cort1sig(cort1count, 1) = arcLength(i,1);
%         cort1sig(cort1count, 3) = i;
%     elseif(drill_states(i,1) == 2)
%         cort2count = cort2count + 1;
%         cort2sig(cort2count, 2) = signature(i,1);
%         cort2sig(cort2count, 1) = arcLength(i,1);
%         cort2sig(cort2count, 3) = i;
%     end
% end
% 
% 
remPts = [];
retPts = [];
for i = 1:1:length(curve(:,1))
    if(signature(i) < 0.5)
        retPts = [retPts;curve(i,:)];
    else
        remPts = [remPts;curve(i,:)];
    end
end


% f1 = figure;
% view(3);
% plot3(curve(:,1), curve(:,3), curve(:,2),'r');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% hold on;
% plot3(smCurve(:,3), smCurve(:,1), smCurve(:,2),'b');
% if(~isempty(remPts))
%     plot3(remPts(:,3), remPts(:,1), remPts(:,2),'ok');
% end
% if(~isempty(retPts))
%     plot3(retPts(:,3), retPts(:,1), retPts(:,2),'og');
% end
% daspect([1,1,1]);
% saveas(f1,'Curve_Kellam_OG','png');
% saveas(f1,'OBCurve_FK_3','png');
% 
% 
% f2 = figure;
% xlabel('Point Indices');
% ylabel('Signature Value');
% plot(cort1sig(:,3),cort1sig(:,2),'Color','#A2142F','LineWidth',2);
% hold on;
% plot(cort2sig(:,3),cort2sig(:,2),'Color','#0072BD','LineWidth',2);
% saveas(f2,'Plot_User_yb','png');
% % saveas(f2,'OBPlot_FK_3','png');
% hold off
% 
% f3 = figure;
% xlabel('Arc Length');
% ylabel('Signature Value');
% plot(cort1sig(:,1), cort1sig(:,2),'Color','#A2142F','LineWidth',2);
% hold on;
% plot(cort2sig(:,1), cort2sig(:,2),'Color','#0072BD','LineWidth',2);
% saveas(f3,'PlotAL_User_yb','png');
% % saveas(f3,'OBPlotAL_FK_3','png');
% hold off
% 
% f4 = figure;
% histogram(signature(:,1));
% saveas(f4,'Hist_User_yb','png');
% saveas(f3,'OBHist_FK_2','png');

% pdf_sig = pdf('Normal',signature);
% f4 = figure;
% plot(signature, pdf_sig,'LineWidth',2);
% saveas(f4,'PDF_User_05','png');