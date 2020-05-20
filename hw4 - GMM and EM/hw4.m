%% Part 1 - Data Generation
mu1 = [-1 -1];
mu2 = [1 1];
mu = [mu1 ; mu2];                   % Means
sigma1 = [.8 0;0 .8];
sigma2 = [.75 -0.2;-0.2 .6];
sigma = zeros(2,2,2);
sigma(:,:,1) = sigma1;
sigma(:,:,2) = sigma2;              % Covariances
p = [0.7 0.3];                      % Mixing proportions p(z=1)=0.7

gmm = gmdistribution(mu,sigma,p);  % gmm is the gaussian mixture model

Y_1000 = random(gmm,1000);    % generates 1000 plots of the GMM

figure          % ploting only gmm
gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
fcontour(gmPDF,[[-4,4] [-4,4]])
grid on
xlabel('x axis')
ylabel('y axis')
title('Gausian Mixture Model')

figure         % ploting scatter and gmm
scatter(Y_1000(:,1),Y_1000(:,2),'.','b');
hold on
gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('1000 scatter plots from GMM')
grid on
xlabel('x axis')
ylabel('y axis')
hold off


% 
% 
% %% Part 2 - K means
% X_50 = random(gmm,50);      % generates 50 plots of the GMM
% K = 20;
% epsilon = 10^-3;            % to stop the iteration
% 
% centroids = zeros(2,2,K+1);
% rng shuffle
% centroids(:,:,1) = -3.5 + (3.5+3.5)*rand(2,2);  % initialized 2 centroids in a square
% %centroids(:,:,1) = random(gmm,2);  % initialized 2 centroids
% 
% figure
% scatter(X_50(:,1),X_50(:,2),'.','g');
% hold on
% scatter(centroids(1,:,1),centroids(2,:,1),'x', 'k');
% grid on
% xlabel('x axis')
% ylabel('y axis')
% title('2-means - centroids initialized')
% legend({'50 random plots','centroids'},'Location','southeast')
% hold off
% 
% for k = 1:K
%     cost = zeros(50,1);         % cost_i = arg(j)min || x_i - mu_j ||
%     clear group1
%     clear group2
%     g1 = 1;
%     g2 = 1;
%     group1 = zeros(1,2);        % saving groups for different centroids
%     group2 = zeros(1,2);        %
%     for i = 1:50                % for each xi - creating the argmin vector
%         C1= norm((X_50(i,:)' - centroids(1,:,k)'),2);
%         C2= norm((X_50(i,:)' - centroids(2,:,k)'),2);
%         if C1 < C2
%             cost(i) = 1;
%             group1(g1,:) = X_50(i,:);
%             g1 = g1 + 1;
%         else
%             cost(i) = 2;
%             group2(g2,:) = X_50(i,:);
%             g2 = g2 + 1;
%         end
%     end
% 
%     % ploting new centroids
%     figure
%     scatter(group1(:,1),group1(:,2),'.','r');
%     hold on
%     scatter(group2(:,1),group2(:,2),'.','b');
%     hold on
%     scatter(centroids(:,1,k),centroids(:,2,k),'x', 'k');
%     grid on
%     xlabel('x axis')
%     ylabel('y axis')
%     title(sprintf('2-means - interation number #%d', k))
%     legend({'centroid 1 group','centroid 2 group','centroids'},'Location','southeast')
%     hold off
% 
%     % Check lecture 6 to see the k-means algorithem
%     upsum1 = zeros(1,2);
%     dosum1 = 0;
%     upsum2 = zeros(1,2);
%     dosum2 = 0;
%     for i = 1:50
%         if cost(i)==1       % argmin is 1
%             upsum1 = upsum1 + 1*X_50(i,:);
%             dosum1 = dosum1 + 1;
%         else                % argmin is 2
%             upsum2 = upsum2 + 1*X_50(i,:);
%             dosum2 = dosum2 + 1;
%         end
%     end
%     centroids(1,:,k+1) = upsum1 ./ dosum1;  % calculating new centroids
%     centroids(2,:,k+1) = upsum2 ./ dosum2;
% 
%     if norm(centroids(:,:,k+1) - centroids(:,:,k)) < epsilon
%         break;  % change is too low
%     end
% end
% 
% % ploting new centroids
% figure
% scatter(group1(:,1),group1(:,2),'.','r');
% hold on
% scatter(group2(:,1),group2(:,2),'.','b');
% hold on
% scatter(centroids(:,1,k+1),centroids(:,2,k+1),'x', 'k');
% grid on
% xlabel('x axis')
% ylabel('y axis')
% title(sprintf('2-means - interation number #%d', k+1))
% legend({'centroid 1 group','centroid 2 group','centroids'},'Location','southeast')
% hold off
% 
% 
% 

%% Part 3 - EM
Data = random(gmm,1000);    % generates 1000 plots of the GMM
K=100;                       % maximum iteration
epsilon = 10^-3;            % to stop the iteration

%making initial guess
Param = struct();
rng shuffle
Param.mu = -3 + (3+3)*rand(1,2,2);

rng shuffle
side = -1 + (1+1)*rand(1);
rng shuffle
diag = abs(side) + 2*rand(1,2);
sig1 = [diag(1,1), side; side, diag(1,2)];

rng shuffle
side = -1 + (1+1)*rand(1);
rng shuffle
diag = abs(side) + 2*rand(1,2);
sig2 = [diag(1,1), side; side, diag(1,2)];

Param.sigma(:,:,1) = sig1;
Param.sigma(:,:,2) = sig2;
rng shuffle
p = rand(1);
Param.phi = [p,1-p];        % probabilities for each gaussian

logloss = [];
for k=1:K
    % Expectation step
    Qz = expectation(Data,Param, 2);  % 2 gaussians. Data_new has the labels
    % Maximization step
    Param = maximization(Qz,Data,Param, 2);% 2 gaussians. Param_new has the new gaussians
    % Log likelihood
    logloss(k) = loglike(Data, Param, 2); % 2 gausians. logloss hopefuly maximized
    
    if k >= 2 && norm(logloss(k)-logloss(k-1))<= epsilon
        break;
    end
end

mu = [Param.mu(:,:,1);Param.mu(:,:,2)];
sigma = Param.sigma;
phi = Param.phi;
gmm = gmdistribution(mu,sigma,phi);  % gmm is the gaussian mixture model

figure          % ploting only gmm
gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
fcontour(gmPDF,[[-4,4] [-4,4]])
grid on
xlabel('x axis')
ylabel('y axis')
title('EM - 2-GMM Only')

figure         % ploting scatter and gmm
scatter(Data(:,1),Data(:,2),'.','b');
hold on
gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('EM - 2-GMM over samples')
xlabel('x axis')
ylabel('y axis')
grid on
hold off

figure          % ploting log likelihood
iter = 1:k;
plot (iter,logloss,'b--o');
title('EM - 2-GMM - Log-liklihood Function')
xlabel('Iteration')
ylabel('Log-Loss')
grid on


%% Part 3 - for 3 
K=100;                       % maximum iteration
epsilon = 10^-3;            % to stop the iteration

%making initial guess
Param = struct();
rng shuffle
Param.mu = -3 + (3+3)*rand(1,2,3);

rng shuffle
side = -1 + (1+1)*rand(1);
rng shuffle
diag = abs(side) + 2*rand(1,2);
sig1 = [diag(1,1), side; side, diag(1,2)];

rng shuffle
side = -1 + (1+1)*rand(1);
rng shuffle
diag = abs(side) + 2*rand(1,2);
sig2 = [diag(1,1), side; side, diag(1,2)];

rng shuffle
side = -1 + (1+1)*rand(1);
rng shuffle
diag = abs(side) + 2*rand(1,2);
sig3 = [diag(1,1), side; side, diag(1,2)];

Param.sigma(:,:,1) = sig1;
Param.sigma(:,:,2) = sig2;
Param.sigma(:,:,3) = sig3;
rng shuffle
Param.phi = [1/3,1/3,1/3];        % probabilities for each gaussian

logloss = [];
for k=1:K
    % Expectation step
    Qz = expectation(Data,Param, 3);  % 3 gaussians. Data_new has the labels
    % Maximization step
    Param = maximization(Qz,Data,Param, 3);% 3 gaussians. Param_new has the new gaussians
    % Log likelihood
    logloss(k) = loglike(Data, Param, 3); % 3 gausians. logloss hopefuly maximized
    
    if k >= 2 && norm(logloss(k)-logloss(k-1))<= epsilon
        break;
    end
end

mu = [Param.mu(:,:,1);Param.mu(:,:,2);Param.mu(:,:,3)];
sigma = Param.sigma;
phi = Param.phi;
gmm = gmdistribution(mu,sigma,phi);  % gmm is the gaussian mixture model

figure          % ploting only gmm
gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
fcontour(gmPDF,[[-4,4] [-4,4]])
grid on
xlabel('x axis')
ylabel('y axis')
title('EM - 3-GMM Only')

figure         % ploting scatter and gmm
scatter(Data(:,1),Data(:,2),'.','b');
hold on
gmPDF = @(x1,x2)reshape(pdf(gmm,[x1(:) x2(:)]),size(x1));
g = gca;
fcontour(gmPDF,[g.XLim g.YLim])
title('EM - 3-GMM over samples')
xlabel('x axis')
ylabel('y axis')
grid on
hold off

figure          % ploting log likelihood
iter = 1:k;
plot (iter,logloss,'b--o');
title('EM - 3-GMM - Log-liklihood Function')
xlabel('Iteration')
ylabel('Log-Loss')
grid on

