%% Construct the input
input = [datamat(1:1442,1:198) label(:)];

[trainedClassifier, validationAccuracy] = trainClassifier(input);
disp("Accuracy: " + validationAccuracy);

%% Try another dataset
%input2 = [datamat(174:1492,1:198)];
yfit = trainedClassifier.predictFcn(trainMat);

actualAnswers = label(trainTS);

correctAnswers = sum(yfit == actualAnswers);
total = length(yfit);
accuracy = correctAnswers/total;

disp("After trained: Accuracy: " + accuracy);

%saveCompactModel(trainedClassifier.ClassificationSVM, 'mySVM')



%% Load the saved SVM
function label = classifyX (X) %#codegen 
%CLASSIFYX Classify using SVM Model 
%  CLASSIFYX classifies the measurements in X 

%  returns class labels in label.

CompactMdl = loadCompactModel('mySVM'); 
label = predict(CompactMdl,X); 
end
