function[mmod,IQRmat,gpred]=posterior(data,opfun,results,mref,Etol)
% POSTERIOR calculates:
% 1) mmod or median of the set of equivalent models (those that
% adjust the data with the same error misfit)
% 2) Inter Quartile Range of the set of equivalent models
% 3) gpred or gravity predicted for each equivalent model
% Authors: Zulima Fernández-Muñiz and Juan Luis Fernández-Martínez
igoodm = find(results.error_hist <=Etol);
good_models = results.historia(igoodm,:);
% solution with good_models
matref = [];
gpred = [];
for k=1:length(igoodm)
    modelk = good_models(k,:);
    [gprek{k},mrefk{k}] = solve_model(modelk,opfun);
    matref = [matref mrefk{k}(:)];
    gpred = [gpred gprek{k}];
end
matref = matref';
gpred = gpred';
medmod = median(matref);
IQRmod = iqr(matref);
IQRmod(IQRmod> prctile(IQRmod,99)) = prctile(IQRmod,99);
IQRmat = reshape(IQRmod,size(mref));
mmod = reshape(medmod, size(mref));
