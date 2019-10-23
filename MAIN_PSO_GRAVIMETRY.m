% MAIN_PSO_GRAVIMETRY solves the forward and inverse synthetic problems, 
% providing graphical representation of:
% 1) median and IQR of the equivalent models with an error misfit less than
% Etol.
% 2) true and best models 
% 3) observed and predicted (from the best model) gravity anomalies
% Authors:  Juan Luis Fernandez-Martinez and Zulima Fernandez-Muniz 2019
% Contact: Juan Luis Fernandez Martinez Email: jlfm@uniovi.es
%===================================================================================
% Forward_Model: This routine contains the observed data.
%===================================================================================
Forward_Model;
%===================================================================================
% Search Space
%===================================================================================
% model.lowlimit   = [rhobmin rhok1min xa11min xa12min za11min za21min ...
%                         rhok2min xa21min xa22min za21min za22min];
% model.upperlimit = [rhobmax rhok1max xa11max xa12max za11max za12max ...
%                         rhok2max xa21max xa22max za21max za22max];
%
model.lowlimit   = [2500 2000 100 200 40  90,  500 400 500  60  85];
model.upperlimit = [3000 4000 250 400 90 180, 2000 800 900 120 190];
%===================================================================================
% PSO options
%===================================================================================
% THE PSO FAMILY
%     case 1 'CP-PSO'   (explorative)
%     case 2 'CC-PSO'   (exploitative)
%     case 3 'PSO'      (exploitative)
%     case 4 'PP-PSO'   (mixte)
%     case 5 'RR_PSO'   (explorative-mixte)
%     case 6 'RC_PSO'   (explorative)
%     case 7 'RP_PSO'   (explorative)
%     case 8 'PR_PSO'   (explorative)
%     case 9 'PC_PSO'   (explorative)
%====================================================================================
algo_type = 5;
PSO_options;
options.pso.maxiter = 50; % maximum number of (input) iterations (only stop criteria)
options.pso.size = 50; % number of particles in each iteration (input)
options.pso.elitism = 0;
%====================================================================================
% OPFUN: problem dependent parameters structure
%====================================================================================
opfun.norm = 2;
opfun.modellog = 0;
opfun.prior.model = []; % Prior model if any
%====================================================================================
% Solving the forward problem
%====================================================================================
funobj = @fcost;
[results] = pso_family(funobj,model,data,options,opfun);
%====================================================================================
% Saving results
%====================================================================================
[i,jmodel] = min(results.error_hist);
if opfun.modellog == 1
    results.historia = 10.^results.historia;
    results.fittest = 10.^results.fittest;
    results.parent = 10.^results.parent;
end
model.inverted = results.historia(jmodel,:)
%====================================================================================
% Solution with the best model (one with less error misfit)
%====================================================================================
[gpre,mref] = solve_model(model.inverted,opfun);
name_results = ['Synthetic',date];
save(name_results)
%====================================================================================
% Posterior analysis
%====================================================================================
Etol = 0.1; % level of misfit error (10% of error misfit)
[mmod,IQRmat,gpred] = posterior(data,opfun,results,mref,Etol);
%====================================================================================
% Graphical representations
%====================================================================================
graphics(rgrid,data,mmod,IQRmat,mref,gpre)