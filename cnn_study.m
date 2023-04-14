clc
clear

%data
trainNoise1 = importdata('bin_train1.mat');
trainTruth1 = importdata('bin_truth1.mat');

trainNoise1 = single(permute(trainNoise1,[1 2 3 5 4]));
trainTruth1 = single(permute(trainTruth1,[1 2 3 5 4]));

% trainNoise2 = importdata('bin_train2.mat');
% trainTruth2 = importdata('bin_truth2.mat');
% 
% trainNoise2 = single(permute(trainNoise2,[1 2 3 5 4]));
% trainTruth2 = single(permute(trainTruth2,[1 2 3 5 4]));
% 
% trainNoise = cat(5,trainNoise1,trainNoise2);
% trainTruth = cat(5,trainTruth1,trainTruth2);
trainNoise = trainNoise1;
trainTruth = trainTruth1;

validateNoise = trainNoise(:,:,:,:,2001:end);
validateTruth = trainTruth(:,:,:,:,2001:end);

trainNoise = trainNoise(:,:,:,:,1:2000);
trainTruth = trainTruth(:,:,:,:,1:2000);

layers = [
    image3dInputLayer([116 9 3],"Name","imageinput")
    convolution3dLayer([9 9 3],32,"Name","conv","Padding","same")
    reluLayer("Name","relu")
    convolution3dLayer([9 9 3],32,"Name","conv_1","Padding","same")
    %batchNormalizationLayer("Name","batchnorm")
    reluLayer("Name","relu_1")
    convolution3dLayer([9 9 3],32,"Name","conv_2","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_1")
    reluLayer("Name","relu_2")
    convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_3")
    convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_3")
    convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_3")
    convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_3")
    convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_3")
    convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
    %batchNormalizationLayer("Name","batchnorm_2")
    reluLayer("Name","relu_3")
%     convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
%     %batchNormalizationLayer("Name","batchnorm_2")
%     reluLayer("Name","relu_3")
%     convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
%     %batchNormalizationLayer("Name","batchnorm_2")
%     reluLayer("Name","relu_3")
%     convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
%     %batchNormalizationLayer("Name","batchnorm_2")
%     reluLayer("Name","relu_3")
%     convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
%     %batchNormalizationLayer("Name","batchnorm_2")
%     reluLayer("Name","relu_3")
%     convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
%     %batchNormalizationLayer("Name","batchnorm_2")
%     reluLayer("Name","relu_3")
%     convolution3dLayer([9 9 3],32,"Name","conv_3","Padding","same")
%     %batchNormalizationLayer("Name","batchnorm_2")
%     reluLayer("Name","relu_3")
    convolution3dLayer([9 9 3],1,"Name","conv_4","Padding","same")
    regressionLayer("Name","regressionoutput")];


miniBatchSize  = 50;
validationFrequency = floor(length(validateTruth)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',20, ...
    'InitialLearnRate',1e-5,...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',50, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{validateNoise,validateTruth}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress', ...
    'Verbose',true);

net = trainNetwork(trainNoise,trainTruth,layers,options);
%%
validatePrediction = predict(net,validateNoise(:,:,:,:,5));
plot(squeeze(validatePrediction(:,5,:)),'o-');
hold on
plot(squeeze(validateNoise(:,5,:,:,5)),'+-');
plot(squeeze(validateTruth(:,5,:,:,5)),'--')






