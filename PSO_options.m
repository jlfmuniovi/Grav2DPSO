% random generation of the initial population (it could be 'given')
options.inversion.seed='random';
options.pso.deltat1=1.2;
options.pso.deltat2=0.8;
options.pso.ccontrol=0; % delta control of coordinates 1: YES
options.pso.cdelta=0.5; % deltat for control
options.pso.proyection='near';   %'far' 'rebond' kind of proyection over search space
options.pso.sinmejora=2;
options.pso.pinterior=50; %percentage
options.inversion.modellog = 0; %1=log or 0=non log parameters for sampling
% deltat for PSO optimization
opfun.deltatmin = 0.1;
opfun.deltatmax = 1.2;
% Control of the time step.
opfun.deltat_option='randp';
% also options.inversion.seed='given'; fro rebooting purposes
opfun.niter=2;
%for lime and sand
opfun.delta1=0.4;
opfun.delta2=0.9;
%==========================================================================
% PSO default parameters
%==========================================================================
% Algorithm type
switch algo_type
    case 1 %%'CP'
        options.pso.esquema = 'CP';
    case 2 %% 'CC'
        options.pso.esquema = 'CC';
    case 3 %%'PSO'
        options.pso.esquema = 'PSO';
    case 4 %%'PP'
        options.pso.esquema = 'PP';
    case 5 %%'RR'
        options.pso.esquema = 'RR';
        % deltat for PSO optimization
        opfun.deltatmin = 0.2;
        opfun.deltatmax = 3;
    case 6 %%'RC'
        options.pso.esquema = 'RC';
    case 7 %%'RP'
        options.pso.esquema = 'RP';
    case 8 %%'PR'
        options.pso.esquema = 'PR';
    case 9 %%'PC'
        options.pso.esquema = 'PC';
    case 10
        options.pso.esquema = 'RN';
end
%--------------------------------------------------------------
% READING THE PARTICLE CLOUDS
%--------------------------------------------------------------
% We swith on the cloud modality
options.pso.cloud='yes';
if options.pso.cloud=='yes'
    if isequal(options.pso.esquema,'PSO')
        disp('loading cloud_PSO.mat');
        load('cloud_PSO.mat');
    elseif isequal(options.pso.esquema,'CP')
        disp('loading cloud_CP.mat');
        load('cloud_CP.mat');
    elseif isequal(options.pso.esquema,'CC')
        disp('loading cloud_CC.mat');
        load('cloud_CC.mat');
    elseif isequal(options.pso.esquema,'PP')
        disp('loading cloud_PP.mat');
        load('cloud_PP.mat');
    elseif isequal(options.pso.esquema,'RR')
        disp('loading cloud_RR.mat');
        load('cloud_RR.mat');
    elseif isequal(options.pso.esquema,'RC')
        disp('loading cloud_RC.mat');
        load('cloud_RC.mat');
    elseif isequal(options.pso.esquema,'RP')
        disp('loading cloud_RP.mat');
        load('cloud_RP.mat');
    elseif isequal(options.pso.esquema,'PR')
        disp('loading cloud_PR.mat');
        load('cloud_PR.mat');
    elseif isequal(options.pso.esquema,'PC') | isequal(options.pso.esquema,'RN') 
        disp('loading cloud_PC.mat');
        load('cloud_PC.mat');
    else
        % for RR the points are calculated analytically
    end
    % loading the cloud
    npoints = size(w_al_ag,1);
    prow = randperm(npoints);
    for i=1:npoints
        dat_w_al_ag(i,1) = w_al_ag(prow(i),1); %w
        dat_w_al_ag(i,2) = w_al_ag(prow(i),2); %al
        dat_w_al_ag(i,3) = w_al_ag(prow(i),3); %ag
    end
    % charging the options structure
    if isequal(options.pso.esquema,'PSO')
        options.pso.pso.inertia=dat_w_al_ag(:,1)';
        options.pso.pso.philocal=dat_w_al_ag(:,2)';
        options.pso.pso.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'CP')
        options.pso.cp.inertia=dat_w_al_ag(:,1)';
        options.pso.cp.philocal=dat_w_al_ag(:,2)';
        options.pso.cp.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'CC')
        options.pso.cc.inertia=dat_w_al_ag(:,1)';
        options.pso.cc.philocal=dat_w_al_ag(:,2)';
        options.pso.cc.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'PP')
        options.pso.pp.inertia=dat_w_al_ag(:,1)';
        options.pso.pp.philocal=dat_w_al_ag(:,2)';
        options.pso.pp.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'RR')
        options.pso.rr.inertia=dat_w_al_ag(:,1)';
        options.pso.rr.philocal=dat_w_al_ag(:,2)';
        options.pso.rr.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'RC')
        options.pso.rc.inertia=dat_w_al_ag(:,1)';
        options.pso.rc.philocal=dat_w_al_ag(:,2)';
        options.pso.rc.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'RP')
        options.pso.rp.inertia=dat_w_al_ag(:,1)';
        options.pso.rp.philocal=dat_w_al_ag(:,2)';
        options.pso.rp.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'PR')
        options.pso.pr.inertia=dat_w_al_ag(:,1)';
        options.pso.pr.philocal=dat_w_al_ag(:,2)';
        options.pso.pr.phiglobal=dat_w_al_ag(:,3)';
    elseif isequal(options.pso.esquema,'PC') | isequal(options.pso.esquema,'RN') 
        options.pso.pc.inertia=dat_w_al_ag(:,1)';
        options.pso.pc.philocal=dat_w_al_ag(:,2)';
        options.pso.pc.phiglobal=dat_w_al_ag(:,3)';
    else
        % no other method
    end
end