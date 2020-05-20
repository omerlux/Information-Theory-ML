% =========================================
%       Homework on K-Nearest Neighbors
% =========================================
% Course: Introduction to Information Theory
% Lecturer: Haim H. Permuter.
%
% NOTE:
% -----
% Please change the variable ID below to your ID number as a string.
% Please do it now and save this file before doing the assignment
clear all;
ID = 'ID';

%% Loading and plot a sample of the data
% ---------------------------------------
% The MNIST database contains in total of 60000 train images and 10000 test images
% of handwritten digits.
% This database is a benchmark for many classification algorithms, including neural networks.
% For further reading, visit http://yann.lecun.com/exdb/mnist/
%
% You will implement the KNN algorithm to classify two digits from the database: 3 and 5.
% First, we will load the data from .MAT file we have prepared for this assignment: MNIST_3_and_5.mat
load('MNIST_3_and_5.mat')
% The data is divided into 2 pairs:
% (Xtrain, Ytrain) , (Xvalid, Yvalid)
% In addition, you have unlabeled test sample in Xtest.
%
% Each row of a X matrix is a sample (gray-scale picture) of dimension 784 = 28^2,
% and the digit (number) is in the corresponding row of the Y vector.
%
% To plot the digit, see attached function 'plot_sample.m'
%%sampleNum = 1;
%%plot_sample(Xvalid(sampleNum,:),Yvalid(sampleNum))


% Build a KNN classifier based on what you've learned in class:
%
% 1. The classifier should be given a train dataset (Xtrain, Ytain),  and a test sample Xtest.
% 2. The classifier should return a class for each row in test sample in a column vector Ytest.
%
% Finally, your results will be saved into a <ID>.txt file, where each row is an element in Ytest.
%
% Note:
% ------
% For you conveniece (and ours), give you classifications in a 1902 x 1 vector named Ytest,
% and set the variable ID at the beginning of this file to your ID.


%% < your code here >
% 
% %% VALID
% error_K_L = zeros(9,4);
% for K=1:2:9  % K neighbors
%     for L=1:4 % L norm       
%         k_neighbors = zeros( K ,2); % a vector of the KNN - {#sample, distance}
%         valid_results = zeros (1522, K+1);   % the first colomn is the valid true value
%         all_strikes = 0;
%         median_strikes =0;
%         
%         %%% VALID
%         for j=1:1522                    % j is for the VALID
%             k_neighbors = zeros( K ,2); % a vector of the KNN - {#sample, distance}
%             d = vecnorm( Xvalid(j,:) - Xtrain , L, 2); %this does a norm on all the Xtrain examples
%             for i=1:11552                % i is for the TRAIN
%                 if i<=K % first K will get in
%                     k_neighbors(i,:) = [d(i), i];
%                 else     % for the rest - replace with the minimum
%                     k_neighbors = sortrows(k_neighbors); % max will be last
%                     if d(i) < k_neighbors(K,1) % if current distance is less than the maximum
%                         k_neighbors(K,:) =[d(i),i]; % replacing for minimum [distance,index of the train]
%                     end
%                 end
%             end
%             
%             %% adding the K neighbors to the array
%             valid_results(j,1) = Yvalid(j); % the true value of the Xvalid is first
%             for k=1:K
%                 index_in_train = k_neighbors (k,2);
%                 valid_results (j,k+1)=Ytrain(index_in_train); % each row has K elements that picked
%             end
%             
%             %%counting strikes
%             m = median(valid_results(j,2:K+1));
%             if (m ~= valid_results(j,1))
%                 median_strikes = median_strikes+1;
%             end
%             for k=2:K
%                 if( valid_results(j,k) ~= valid_results(j,1))
%                     all_strikes = all_strikes+1;
%                 end
%             end
%         end
%         error_K_L(K,L)= median_strikes; % strikes from 1522 valid tests
%     end
% end
% %% for L4,K=5 17/1522 are wrong

%%% TEST
%% creating the Ytest vector
K=5;
L=4;

Ytest = zeros(1902,1);
test_results = zeros(1902,K);

k_neighbors = zeros( K ,2); % a vector of the KNN - {#sample, distance}
all_strikes = 0;
median_strikes =0;

for j=1:1902                    % j is for the TEST
    k_neighbors = zeros( K ,2); % a vector of the KNN - {#sample, distance}
    d = vecnorm( Xtest(j,:) - Xtrain , L, 2); %this does a norm on all the Xtrain examples
    for i=1:11552                % i is for the TRAIN
        if i<=K % first K will get in
            k_neighbors(i,:) = [d(i), i];
        else     % for the rest - replace with the minimum
            k_neighbors = sortrows(k_neighbors); % max will be last
            if d(i) < k_neighbors(K,1) % if current distance is less than the maximum
                k_neighbors(K,:) =[d(i),i]; % replacing for minimum [distance,index of the train]
            end
        end
    end
    
    %% adding the K neighbors to the array
    for k=1:K
        index_in_train = k_neighbors (k,2);
        test_results (j,k)=Ytrain(index_in_train); % each row has K elements that picked
    end
    
    %% creating the Ytest results
    m = median(test_results(j,:));
    Ytest(j) = m;
end
   

%% save classification results
disp('saving')
csvwrite([ID '.txt'], Ytest)
disp('done')

% %% Trying to print results
% index = randperm(1902);
% for j=1:1902
%     cur_index = index(j);
%     plot_sample(Xtest(cur_index,:),Ytest(cur_index))
% end