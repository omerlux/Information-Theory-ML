function LogLike = loglike(Data, Param, numG)
%returning the log likelihood of the function - eq.10b
LogLike = 0;

if numG == 3
    for i = 1:size(Data,1)
        x = Data(i,1:2);  % taking 1 sample
        % gi is the prob of gaussian x sample
        g1 = Param.phi(1,1) * mvnpdf(x, Param.mu(:,:,1), Param.sigma(:,:,1));
        g2 = Param.phi(1,2) * mvnpdf(x, Param.mu(:,:,2), Param.sigma(:,:,2));
        g3 = Param.phi(1,3) * mvnpdf(x, Param.mu(:,:,3), Param.sigma(:,:,3));
        LogLike = LogLike + log(g1+g2+g3);
    end
else
    for i = 1:size(Data,1)
        x = Data(i,1:2);  % taking 1 sample
        % gi is the prob of gaussian x sample
        g1 = Param.phi(1,1) * mvnpdf(x, Param.mu(:,:,1), Param.sigma(:,:,1));
        g2 = Param.phi(1,2) * mvnpdf(x, Param.mu(:,:,2), Param.sigma(:,:,2));
        LogLike = LogLike + log(g1+g2);
    end
end
end
