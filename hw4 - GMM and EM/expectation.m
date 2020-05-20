function Qz = expectation(Data, Param, numG)
% returns Qz - the weights Q which the gaussians will be computed by

mu1 = Param.mu(:,:,1);
mu2 = Param.mu(:,:,2);
sigma1 = Param.sigma(:,:,1);
sigma2 = Param.sigma(:,:,2);
phi = Param.phi;

if numG == 3
    mu3 = Param.mu(:,:,3);
    sigma3 = Param.sigma(:,:,3);
    Qz = zeros(size(Data,1),3);
else
    Qz = zeros(size(Data,1),2);
end


for i=1:size(Data,1)
    x = Data(i,1:2);  % taking 1 sample
    
    % mvnpdf will be the W_ij and the phi is the phi(j)
    p1 = phi(1,1) * mvnpdf (x, mu1, sigma1);   % check probability of gaussian 1
    p2 = phi(1,2) * mvnpdf (x, mu2, sigma2);   % check probability of gaussian 2
    
    % Qz is according to 29a,b
    if numG == 3
        p3 = phi(1,3) * mvnpdf (x, mu3, sigma3);   % check probability of gaussian 3
        Qz(i,1:3) = [p1/(p1+p2+p3), p2/(p1+p2+p3), p3/(p1+p2+p3)]; % creating Q - weights
    else
        Qz(i,1:2) = [p1/(p1+p2), p2/(p1+p2)];      % creating Q - weights
    end
end

