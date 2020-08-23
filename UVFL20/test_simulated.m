%%************************************************************************
%% run random matrix completion problems. 
%% ************************************************************************

clear all;

restoredefaultpath;

addpath(genpath('solvers'));

addpath(genpath('PROPACKmod'));
%%

OPTIONS_AMM.maxiter = 5000;

OPTIONS_AMM.printyes = 1;

OPTIONS_MAPM.maxiter = 2000;

OPTIONS_MAPM.printyes = 1;

OPTIONS_ALS.maxiter = 2000;

OPTIONS_ALS.printyes = 1;

OPTIONS_ALS.tol = 1.0e-6;

%% generate random a test problem

scenario = 'noisy';

ntest = 5;       

%% ************** Initialization for test problem ******************

nr = 1000;  nc = 1000;     % nr=3000,  nr=5000

rstar = 10;   % 10   20

SR = [0.1   0.15    0.2   0.25];

ns = length(SR);

nuAMM  = [45    40    40     40]; 

nuMAPM = [10    10    10     10];

nuALS = [1    2.5    1.8    1.5];

%% ***************** Initialization *********************************

AMM_matrelerr = zeros(ntest,ns); Brid_matrelerr = zeros(ntest,ns);  ALS_matrelerr = zeros(ntest,ns);

AMM_matrank = zeros(ntest,ns);   Brid_matrank = zeros(ntest,ns);    ALS_matrank = zeros(ntest,ns);

AMM_mattime = zeros(ntest,ns);   Brid_mattime = zeros(ntest,ns);    ALS_mattime = zeros(ntest,ns);

%% ******************** main loop  **********************************************

for i = 1:4
    
        const_AMM  = nuAMM(i);
    
        const_MAPM = nuMAPM(i);
    
        const_ALS  = nuALS(i);
    
    for test_iter = 1:ntest
        
        test_iter
        
        randstate =  100*i*test_iter
        randn('state',double(randstate));
        rand('state',double(randstate));
        
        if strcmp(scenario,'noiseless')
            noiseratio = 0;
        else
            noiseratio = 0.1;
        end
        
        fprintf('\n nr = %2.0d,   nc = %2.0d,   rank = %2.0d\n',nr,nc,rstar);
        randn('state',double(randstate));
        rand('state',double(randstate));
        
        p = round(SR(i)*nr*nc);     %% number of sampled entries
        
        %% *************** to generate the true matrix ******************
        
        M.U = randn(nr,rstar);
        
        M.V = randn(nc,rstar);
        
        normM = sqrt(sum(sum((M.U'*M.U).*(M.V'*M.V))));
        
        Mstar = M.U*M.V';
        
        num_sample = p;
        
    %%  *************  non-uniform sampling  **********************
        fprintf('\n non-uniform sampling\n');
        fprintf('\n SR = %2.2f, noiseratio = %2.2f\n,',SR(i),noiseratio);
        pvec = ones(nr,1);
        cnt = round(0.1*nr);
        pvec(1:cnt) = 2*pvec(1:cnt);
        pvec(cnt+[1:cnt]) = 4*pvec(cnt+[1:cnt]);
        pvec = nr*pvec/sum(pvec);
        qvec = ones(nc,1);
        cnt = round(0.1*nc);
        qvec(1:cnt) = 2*qvec(1:cnt);
        qvec(cnt+[1:cnt]) = 4*qvec(cnt+[1:cnt]);
        qvec = nc*qvec/sum(qvec);
        probmatrix = rand(nr,nc).*(pvec*qvec');
        [probvec,sortidx] = sort(probmatrix(:),'descend');
        nzsortidx = find(probvec>= probvec(p));
        nzidx = sortidx(nzsortidx);
        zidx = sortidx(p+1:end);
             
        bb =  Mstar(nzidx);
        
        if strcmp(scenario,'noiseless')
            xi = sparse(p,1);
            sigma = 0;
        else
            randnvec = randn(p,1);
            sigma = noiseratio*norm(bb)/norm(randnvec);
            xi = sigma*randnvec;
            bb = bb + xi;
        end
        
        A = zeros(nr,nc);
        
        A(nzidx) = bb;
        
     %% ***************** Initialization part ************************
        
        mu = 1.0e-8;
        
        normb = norm(bb);
        
        r = min(min(nr,nc),150);
        
        pars.normM = normM;
        
        pars.normb = normb;
        
        pars.nc = nc;  pars.nr = nr;
        
       %% ********************** to seek the starting point ****************************** 
        tstart = clock; 
       
        options = 1.0e-6;
     
        [U,dA,V] = lansvd(full(A),r,'L',options);
       
        Ustart = U(:,1:r);
        
        Vstart = V(:,1:r);
        
        time01 = etime(clock,tstart);
   
     %% **************** used for MAPM and ALS *********************
        dd = ones(1,r);
        
        UVstart = (Ustart.*dd)*Vstart';
        
        Mstart = zeros(nr,nc);
        
        Mstart(nzidx) = A(nzidx);
        
        Mstart(zidx) = UVstart(zidx);
        
        time02 = etime(clock,tstart);
        
        dA = diag(dA)';
        
        max_dA = max(dA);
        
  %% ********************** AMM_solver ******************************
        
        OPTIONS_AMM.Lip_const = 2.5*max_dA;
        
        lambda = 10*const_AMM*normb*SR(i);

        OPTIONS_AMM.tol = 1.0e-3;
        
        tstart = clock;
        
        Ustart_AMM = Ustart.*dA(1:r).^(1/2);  %
        
        Vstart_AMM = Vstart.*dA(1:r).^(1/2);  %
        
        clear U  V  dA;
        
        [AMM_Xopt,AMM_rank] = AMM_solver(A,Ustart_AMM,Vstart_AMM,nzidx,OPTIONS_AMM,pars,lambda,mu,r);
        
        AMM_time = etime(clock,tstart) + time01;
        
        AMM_relerr = norm(AMM_Xopt- Mstar,'fro')/normM;
        
        AMM_matrank(test_iter,i) = AMM_rank;
        
        AMM_matrelerr(test_iter,i) = AMM_relerr;
        
        AMM_mattime(test_iter,i) =  AMM_time;
         
        %% ********************** Hybrid_solver *************************
        
        OPTIONS_AMM.Lip_const = 2.5*max_dA;
        
        lambda = 10*const_MAPM*normb*SR(i);
        
        tstart = clock;
        
        [Uinit,Vinit,Brid_rank]= Hybrid_initial(Mstart,Ustart,Vstart,dd,zidx,nzidx,OPTIONS_MAPM,pars,lambda,mu,r);
        
        OPTIONS_AMM.tol = 5.0e-3;
        
        Xopt = Hybrid_smooth(A,Uinit,Vinit,nzidx,OPTIONS_AMM,pars,mu);
        
        Brid_time = etime(clock,tstart) + time02;
        
        Brid_relerr = norm(Xopt- Mstar,'fro')/normM;
        
        Brid_matrelerr(test_iter,i) = Brid_relerr;
        
        Brid_matrank(test_iter,i) = Brid_rank;
        
        Brid_mattime(test_iter,i) = Brid_time;
        
     %% ********************** ALS_solver ******************************
        
        lambda = const_ALS*max_dA*SR(i);
        
        tstart = clock;
        
        [FXopt,Fro_rank] = ALS_solver(Mstart,Ustart,Vstart,dd,zidx,nzidx,OPTIONS_ALS,lambda,r);
        
        Fro_time = etime(clock,tstart) + time02;
        
        Fro_relerr = norm(FXopt- Mstar,'fro')/normM;
        
        ALS_matrelerr(test_iter,i) = Fro_relerr;
        
        ALS_matrank(test_iter,i) = Fro_rank;
        
        ALS_mattime(test_iter,i) = Fro_time;
         
    end
    
    AMM_averelerr(i)= mean(AMM_matrelerr(:,i))
    
    AMM_averank(i)= mean(AMM_matrank(:,i))
    
    AMM_avetime(i)= mean(AMM_mattime(:,i))
    
    
    Brid_averelerr(i) =  mean(Brid_matrelerr(:,i))
    
    Brid_averank(i) =  mean(Brid_matrank(:,i))
    
    Brid_avetime(i) =  mean(Brid_mattime(:,i))
    
    
    ALS_averelerr(i) =  mean(ALS_matrelerr(:,i))
    
    ALS_averank(i) = mean(ALS_matrank(:,i))
    
    ALS_avetime(i) =  mean(ALS_mattime(:,i))
    
end

%% *************************************************************************
