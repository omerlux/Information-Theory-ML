function Param = maximization(Qz, Data, Param, numG)
% Maximization step -
%returns the new parameters of the gaussians according to the labels
% we gathered from the expectation step

% Calculating phi - eq. 30
Param.phi = sum(Qz(:,:)) / size(Qz,1);

% Calculating mu - eq. 31
% mu1
Param.mu(:,:,1) = [sum(Qz(:,1) .* Data(:,1)), sum(Qz(:,1) .* Data(:,2)),] ./ sum(Qz(:,1));
% mu2
Param.mu(:,:,2) = [sum(Qz(:,2) .* Data(:,1)), sum(Qz(:,2) .* Data(:,2)),] ./ sum(Qz(:,2));
% mu3
if numG == 3
    Param.mu(1,:,3) = [sum(Qz(:,3) .* Data(:,1)), sum(Qz(:,3) .* Data(:,2)),] ./ sum(Qz(:,3));
end

% Calculating sigma - eq. 32
% sigma 1
sum_sigma1 = zeros(2,2);
for i=1:size(Data,1)
    sum_sigma1 = sum_sigma1 +  Qz(i,1).* ...
        ((Data(i,:)' - Param.mu(:,:,1)')*(Data(i,:)' - Param.mu(:,:,1)')');
end
Param.sigma(:,:,1) = sum_sigma1 ./ sum(Qz(:,1));
% sigma 2
sum_sigma2 = zeros(2,2);
for i=1:size(Data,1)
    sum_sigma2 = sum_sigma2 +  Qz(i,2).* ...
        ((Data(i,:)' - Param.mu(:,:,2)')*(Data(i,:)' - Param.mu(:,:,2)')');
end
Param.sigma(:,:,2) = sum_sigma2 ./ sum(Qz(:,2));
% sigma 3
if numG == 3
    sum_sigma3 = zeros(2,2);
    for i=1:size(Data,1)
        sum_sigma3 = sum_sigma3 +  Qz(i,3).* ...
            ((Data(i,:)' - Param.mu(:,:,3)')*(Data(i,:)' - Param.mu(:,:,3)')');
    end
    Param.sigma(:,:,3) = sum_sigma3 ./ sum(Qz(:,3));
end
end