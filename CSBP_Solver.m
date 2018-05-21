function [X, varargout] = CSBP_Solver(Y, G, opt)

% To get this to work, we're going to have to paste in all the code from
% these functions because they rely on the fact that MATLAB used to let you
% read in function specific variables like vargin from the global scope if
% called from within a function.  This is no longer the case.
% check_size();
% Simple checks on the sizes of Y and G and varG for bloc matrices

if ((strcmp(opt.method,'AMPseededHadamard') == 0) && (strcmp(opt.method,'AMPseededHadamardTranspose') == 0) && (strcmp(opt.method,'AMPseededHadamardTransposeA') == 0) && (strcmp(opt.method,'AMPseededFourier') == 0) )
    [M,N] = size(G);
    if ((M > N) && (opt.alphaBig ~= 1) )
        G = G.';
        [M,N] = size(G);
    end
else
    M = opt.M; N = opt.N;
end

[a,b] = size(Y);
if (a > b)
    Y = Y.';
end
measure_rate = M ./ N;

if (strcmp(opt.method,'AMPhb') )
    [MM,NN] = size(opt.varG);
    if ((MM ~= M) || (NN ~= N) )
        if (((sum(opt.Nvec_bm) == []) > 0) | ((sum(opt.Mvec_bm) == []) > 0) ); error('varG has not the same size as G and thus opt.Mvec_bm and opt.Nvec_bm must be given as input');
        else opt.varG = convert2fullMN(opt.varG,opt.Mvec_bm,opt.Nvec_bm); end;
    end
    if (MM > NN)
        opt.varG = opt.varG.';
        [MM,NN] = size(opt.varG);
    end
end


% set_parameters();
% Reading parameters
if (nargin <= 2); opt = CSBP_Solver_Opt(); end;

% Set parameters
Gmean = opt.Gmean; Ymean = opt.Ymean;
Nvec = opt.Nvec; if (Nvec == 1); Nvec = N; end;
Mvec = opt.Mvec; if (Mvec == 1); Mvec = M; end;

% Prior dependent quantities
switch opt.prior
    case 'SparseGauss'
        prior_param{1} = opt.m_gauss; prior_param{2} = opt.var_gauss; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SparseGaussCut'
        prior_param{1} = opt.m_gaussC; prior_param{2} = opt.var_gaussC; prior_param{3} = opt.Cut; prior_param{4} = NaN;
    case 'SparseGaussPositive'
        prior_param{1} = opt.m_gaussP; prior_param{2} = opt.var_gaussP; prior_param{3} = NaN; prior_param{4} = NaN;
    case '2Gauss'
        prior_param{1} = opt.m_1_gauss; prior_param{2} = opt.m_2_gauss; prior_param{3} = opt.var_1_gauss; prior_param{4} = opt.var_2_gauss;
    case 'SparseExponential'
        prior_param{1} = opt.expo; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SparseConstant'
        prior_param{1} = opt.c_down; prior_param{2} = opt.c_up; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SparseBinary'
        prior_param{1} = NaN; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'L1'
        prior_param{1} = opt.min; prior_param{2} = opt.max; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'Laplace'
        prior_param{1} = opt.beta; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'Binary1'
        prior_param{1} = NaN; prior_param{2} = NaN; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'Complex'
        prior_param{1} = opt.m_gaussComp; prior_param{2} = opt.var_gaussComp; prior_param{3} = NaN; prior_param{4} = NaN;
    case 'SuperpositionCode'
        prior_param{1} = opt.NbSC; prior_param{2} = opt.nonZeroValues; prior_param{3} = NaN; prior_param{4} = NaN;
    otherwise
        disp('unknown prior')
end

if (opt.signal_rho < 0); opt.signal_rho = measure_rate ./ 10; end;
[a,b] = size(opt.signal); if (a > b); opt.signal = opt.signal.'; end;

if (strcmp(opt.method,'AMP') || strcmp(opt.method,'AMPh') )
    if (opt.save_speed == 1) % requires approximately 3 to 4 times the memory space needed to store the measurement matrix G but speed up the algorithm.
        if (opt.remove_mean == 0); varG = sum(sum(G.^2) ) ./ (M .* N); end;
        if (opt.remove_mean > 0)
            % Remove mean in G and Y
            if (max(size(Gmean) ) < 1)
                % Not provided by user, compute it!
                disp('Computing Mean');
                Gmean = sum(G) ./ M; % [1,N]
                Gmean1d = Gmean; % [1,N]
                if (opt.remove_mean == 2); Gmean = sum(Gmean) ./ N; else Gmean = ones(M,1) * Gmean; end;
            else
                disp('Using provided Means');
            end
            Yeff = Y - sum(Y) ./ M;
            Gm2 = Gmean.^2;
            Gm1d2 = Gmean1d.^2;
            GGm = G .* Gmean;
            Geff = G - Gmean;
        end
        if (strcmp(opt.method,'AMP') ); G2 = G.^2;
        elseif (strcmp(opt.method,'AMPh') && (opt.remove_mean == 2) ); G2m = sum(sum(G.^2) ) ./ (M .* N);
        elseif (strcmp(opt.method,'AMPh') && (opt.remove_mean == 1) ); G2m = sum(G.^2) ./ M;
        end;
    end
    
    if (opt.save_memory == 1) % requires 3 to 4 times less memory but the algorithm is slower because pre-computation are not allowed.
        if (opt.remove_mean == 0); varG = sum(sum(G.^2) ) ./ (M .* N); end;
        if (opt.remove_mean > 0)
            % Remove mean in G and Y
            if (max(size(Gmean) ) < 1)
                % Not provided by user, compute it!
                disp('Computing Mean');
                Gmean = sum(G) ./ M; % [1,N]
                Gmean1d = Gmean; % [1,N]
                if (opt.remove_mean == 2); Gmean = sum(Gmean) ./ N; end;
            else
                disp('Using provided Means');
            end
            Yeff = Y - sum(Y) ./ M;
            Gm2 = Gmean.^2;
            Gm1d2 = Gmean1d.^2;
        end
        if (strcmp(opt.method,'AMPh') && (opt.remove_mean == 2) ); G2m = sum(sum(G.^2) ) ./ (M .* N); varG = abs(G2m - Gm2);
        elseif (strcmp(opt.method,'AMPh') && (opt.remove_mean == 1) ); G2m = sum(G.^2) ./ M; varG = abs(G2m - Gm2);
        end;
    end
end

if (strcmp(opt.method,'AMPhadamard') == 1)
    load Index_random_part1;
    M = opt.M; N = opt.N;
    if (opt.remove_mean > 0); Gmean = 1 ./ M .* Hadamard_multiply_vec_prime(ones(1,M),M,N,Index_random_part1); end
    Gmean1d = Gmean; % [1,N]
    if (opt.remove_mean == 2); Gmean = sum(Gmean) ./ N; end
    Yeff = Y - sum(Y) ./ M;
    Gm2 = Gmean.^2;
    Gm1d2 = Gmean1d.^2;
end

if ((strcmp(opt.method,'AMPseededHadamard') == 1 ) || (strcmp(opt.method,'AMPseededHadamardTranspose') == 1) )
    M = opt.M; N = opt.N;
    if (opt.remove_mean > 0);
        
        vec = zeros(1, opt.M);
        Gmean = zeros(opt.numBlockL, opt.N);
        Ymean = zeros(1, opt.M);
        for l = 1 : opt.numBlockL
            vec(sum(opt.Mblock(1 : (l - 1) ) ) + 1 : sum(opt.Mblock(1 : l) ) ) = 1;
            Gmean(l, :) = 1 / opt.Mblock(l) * MultSeededHadamardTranspose2(vec, opt.J, opt.numBlockL, opt.numBlockC, opt.Mblock, opt.Nblock, opt.rp, opt.noBlockError).';
            Ymean(sum(opt.Mblock(1 : (l - 1) ) ) + 1 : sum(opt.Mblock(1 : l) ) ) = mean(Y(sum(opt.Mblock(1 : (l - 1) ) ) + 1 : sum(opt.Mblock(1 : l) ) ) );
            vec = zeros(1, opt.M);
        end
        
        Gmean2 = Gmean.^2;
        Yeff = Y - Ymean;
    end
end

display_information();
initialisation();

% Construction of the prior dependent class and noise_and_error
prior = Prior(opt.signal_rho, N, measure_rate, opt.learn, opt.prior, opt.dump_learn, R_init, S2_init, av_mess_init, var_mess_init, opt.method, prior_param{1}, prior_param{2}, prior_param{3}, prior_param{4}); F = str2func(prior.func);
n_and_e = noise_and_error(opt.conv, opt.var_noise, opt.dump_learn);

t = 1;
% print infos to screen
if (opt.print >= 1); print_to_screen(); end
% initial error (should be rho)
if (max(size(opt.signal) ) > 2);
    if(strcmp(opt.prior, 'SuperpositionCode') == 0); n_and_e = n_and_e.compute_true_MSE(opt.signal,prior.av_mess);
    else n_and_e = n_and_e.compute_true_SER(opt.signal,prior.av_mess,opt.NbSC,opt.nonZeroValues); end;
    MSEt(t) = n_and_e.true_error;
end
% initial MSE by block
if ((opt.MSEbyBlock > 0) && (mod(t, opt.MSEbyBlock) == 0) && (opt.numBlockC > 1) ); MSEblocks{t} = MSEbyBloc(prior.av_mess, opt.signal, opt.numBlockC, opt.Nblock); end

figure;
% Starting main code
while (t <= opt.nb_iter)
    
    switch (opt.method)
        case ('AMP'); AMP();
        case ('AMPcomplex'); AMPcomplex();
        case ('AMPtap'); AMPtap();
        case ('AMPseededHadamard'); AMPseededHadamard();
        case ('AMPseededFourier'); AMPseededFourier();
    end
    
    % Test of the convergence
    n_and_e = n_and_e.compute_convergence(prior.av_mess_old,prior.av_mess);
    if ((n_and_e.convergence < opt.conv) && (t > 10) ); fprintf('Converged, convergence = %e',n_and_e.convergence); break; end;
    
    % Test of reconstruction on the fly knowing the original signal
    if (max(size(opt.signal) ) > 2)
        if(strcmp(opt.prior, 'SuperpositionCode') == 0); n_and_e = n_and_e.compute_true_MSE(opt.signal,prior.av_mess);
        else n_and_e = n_and_e.compute_true_SER(opt.signal,prior.av_mess,opt.NbSC,opt.nonZeroValues); end;
        if ((n_and_e.true_error < opt.conv) ); fprintf('Solution found, true error = %e\n',n_and_e.true_error); break; end;
    end
    
    % Learning of the noise if activated
    if (opt.option_noise == 1); n_and_e = n_and_e.learn_noise(Y,W_new,V_new); end
    
    % print infos to screen
    if ((opt.print > 0) && (mod(t, opt.print) == 0) ); print_to_screen(); end
    
    % MSE by block
    if ((opt.MSEbyBlock > 0) && (mod(t, opt.MSEbyBlock) == 0) && (opt.numBlockC > 1) )
        MSEblocks{t + 1} = MSEbyBloc(prior.av_mess, opt.signal, opt.numBlockC, opt.Nblock);
        semilogy(MSEblocks{t + 1} );
        if (t == opt.MSEbyBlock); MSEmax = max(max(MSEblocks{t + 1} ) ); end
        axis ([1, opt.numBlockC, 0, MSEmax] ); drawnow;
    end
    
    % MSE as a function of the iterations, to be compared with density evolution
    if (max(size(opt.signal) ) > 2)
        MSEt(t + 1) = n_and_e.true_error;
        if ((isnan(n_and_e.true_error) == 1) || (n_and_e.true_error > 1e5) ); error('The algorithm did not converge'); end
    end
    
    t = t + 1;
    
end

close(gcf);

% can output the MSE a function of iterations
varargout{1} = prior;
varargout{2} = n_and_e;
if (max(size(opt.signal) ) > 2); varargout{3} = MSEt; end
if ((opt.MSEbyBlock > 0) && (opt.numBlockC > 1) ); varargout{4} = MSEblocks; else varargout{4} = []; end

X = prior.av_mess;

end