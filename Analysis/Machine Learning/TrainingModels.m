clc;
clear;

% reading the data
data = xlsread('Data.xlsx');

% shuffling the rows
[m, n] = size(data);
idx = randperm(m);
shuffled_data = data;
shuffled_data(idx,:) = data(:,:);

% initiating the num of features and trees
feature_cols = [1 3 5];
num_features = length(feature_cols);
num_trees = 100;


% initiating the feature and column matrix
features = zeros(length(shuffled_data(:,1)),num_features);
classes = zeros(length(shuffled_data(:,1)),1);

% assigning values to the feature and column matrix
features = shuffled_data(:,feature_cols);
classes = shuffled_data(:,7);

% partitioning data
split = 0.70;
[trainX, trainY, testX, testY] = PartitionData(features, classes, split);

% training data
% training on Random Forests
randomForest=TreeBagger(num_trees,trainX,trainY,'OOBPred','On','OOBVarImp',...
    'on','Method','classification');

% training on k-nn
knn = fitcknn(trainX,trainY,'NumNeighbors',1,'Standardize',0);

% training on Naive Bayes Classifier
nb = fitcnb(trainX, trainY);

% testing data
% Random Forest
[classificationAccuracy_RF, classificationProbab_RF, predictions_RF,predMetrics_RF]...
    = ComputeClassificationAccuracyRF(randomForest,testX,testY);

% K-Nearest Neighbors
[accuracy_knn, predictions_knn] = ComputeClassificationAccuracyGeneral(knn, testX, testY);

% Naive Bayes
[accuracy_cnb, predictions_cnb] = ComputeClassificationAccuracyGeneral(nb, testX, testY);