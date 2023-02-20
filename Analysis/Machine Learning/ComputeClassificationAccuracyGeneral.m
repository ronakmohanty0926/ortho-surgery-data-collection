function [accuracy, predictions] = ComputeClassificationAccuracyGeneral(Model,testX, testY)

groundTruth = testY;

%Predict using the random forest model
testDataNew = testX;
predictions = predict(Model,testX);

%Calculate accuracy
numSuccess = length(find(predictions-groundTruth == 0));
accuracy = numSuccess/length(groundTruth);
