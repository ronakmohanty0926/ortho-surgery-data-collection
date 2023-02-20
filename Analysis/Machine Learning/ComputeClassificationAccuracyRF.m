%This function predicts a given curve using the randomForest model provided during the function call
%Following inputs are needed:
%1. randomForest: This is the trained random forest model which is to be used for predictions
%2. testX:		  specify the curve feature set to be predicted using the trained model. Ensure the format of
%				  testX and the data random forest model is trained on is the same
%3. testY:		  Specify the actual class of testX vectors to compute the model's classification accuracy
%The function returns the prediction accuracy of testX, probabilities of each prediction, 
%the 0/1 predictions, and the prediction metrics (precision, recal, true positive rate, true negative rate)

function [accuracy, classificationProbab, predictions, predictionMetrics] = ComputeClassificationAccuracyRF(randomForest,testX, testY)

groundTruth = testY;

%Predict using the random forest model
testDataNew = testX;
[classification, classificationProbab]= randomForest.predict(testDataNew);

for i=1:length(classification)
    if(classification{i}=='0')
        predictions(i,1) = 0;
    elseif(classification{i}=='1')
        predictions(i,1) = 1;
    elseif(classification{i}=='2')
        predictions(i,1) = 2;
    elseif(classification{i}=='3')
        predictions(i,1) = 3;
    end
end

numSuccess = length(find(predictions-groundTruth == 0));
accuracy = numSuccess/length(groundTruth);

%Identify false positives & false negatives
falsePositives=[]; falseNegatives=[];
p=1; n=1;
for i=1:length(predictions)
    if(groundTruth(i,end)==1 && predictions(i,end)==0)
        falseNegatives(n,1)=i; n=n+1;
    end
    if(groundTruth(i,end)==0 && predictions(i,end)==1)
        falsePositives(p,1)=i; p=p+1;
    end
end

%Identify true positives & true negatives
truePositives=[]; trueNegatives=[];
p=1; n=1;
for i=1:length(predictions)
    if(groundTruth(i,end)==0 && predictions(i,end)==0)
        trueNegatives(n,1)=i; n=n+1;
    end
    if(groundTruth(i,end)==1 && predictions(i,end)==1)
        truePositives(p,1)=i; p=p+1;
    end
end

%Compute precision & recall
precision = length(truePositives)/(length(truePositives) + length(falsePositives));
recall = length(truePositives)/(length(truePositives) + length(falseNegatives));

%True negative & true positive rate
true_neg_rate = length(trueNegatives)/(length(trueNegatives) + length(falsePositives));
true_pos_rate = length(truePositives)/(length(truePositives) + length(falseNegatives));

predictionMetrics = [accuracy precision recall true_neg_rate true_pos_rate];