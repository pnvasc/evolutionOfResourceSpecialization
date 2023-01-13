% How does joint evolution of consumer traits affect resource specialization?
% Paula Vasconcelos & Claus Rueffler
% Matlab code accompanying the paper with the same title published in The American Naturalist in 2019

function EH_parameterSweep_Dryad

clc
rng(2)
rng('shuffle')

%%% INITIALIZE SIMULATION PARAMETERS %%%

iter  = 1;    % Number of simulation iterations
nt    = 1e6;  % Number of consumer generations
h     = 1e-4; % Step size for Euler's method in resource dynamics
Data  = [];   % Empty array where information about evolutionary outcome and parameter combinations that lead to them are collected

%%% INITIALIZE POPULATION PARAMETERS %%%

te    = [];   % Trade-off for feeding efficiency
th    = [];   % Trade-off for handling time
r     = 1e4;
K     = 2.5;
a     = 1;
d     = 0.9;
u     = 0.005; % Mutation rate

%%% MUTATION MATRIX %%%

sige  =  0.01;
sigh  =  0.01;
rho   =  0;
sigma =  [sige^2 sige*sigh*rho; sige*sigh*rho sigh^2];

for     te = -5:0.2:5
    for th = -5:0.2:5
        V  = [te;th];
        
        if  te == 0
            te =  1e-3;
        end
        
        if  th == 0
            th =  1e-3;
        end
        
        for z = 1:iter
            
            %%% INITIALIZE POPULATIONS %%%
            
            R1 = K;
            R2 = K;
            C  = 100;
            
            %%% INITIALIZE TRAITS %%%
            
            xe = 0.35;
            xh = 0.35;
            e1 = [1     * ((1 - exp(1) .^ (-te .* (1 - xe))) ./ (1 - exp(1) .^ (-te)))];
            e2 = [1     * ((1 - exp(1) .^ (-te .*      xe))  ./ (1 - exp(1) .^ (-te)))];
            h1 = [1 - 1 * ((1 - exp(1) .^ (-th .* (1 - xh))) ./ (1 - exp(1) .^ (-th)))];
            h2 = [1 - 1 * ((1 - exp(1) .^ (-th .*      xh))  ./ (1 - exp(1) .^ (-th)))];
            
            
            %%%% SIMULATION %%%%
            
            rng('shuffle')
            for t = 1:nt
                
                %%% RRESOURCE DYNAMICS %%%
                
                start     = 0;
                counter   = 1e1;
                while start < counter
                    R1    = R1 + (r * R1 * (1 - R1 / K) - sum(((e1 .* R1 .*C) ./ (1 + e1 .* h1 .* R1 + e2 .* h2 .* R2)))) * h;
                    R2    = R2 + (r * R2 * (1 - R2 / K) - sum(((e2 .* R2 .*C) ./ (1 + e1 .* h1 .* R1 + e2 .* h2 .* R2)))) * h;
                    start = start + 1;
                end
                
                %%% CONSUMER DYNAMICS %%%
                
                fC = (poissrnd((C .* a .* e1 .* R1 + C .* a .* e2 .* R2) ./ (1 + e1 .* h1 .* R1 + e2 .* h2 .* R2))) - binornd(C, d);
                C  = C + fC;
                
                if isempty(C)
                    continue
                end
                
                if isnan(C)
                    continue
                end
                
                %%% MUTATIONS %%%
                
                mut_off      = binornd(C,u);              % Number of mutants
                mi           = find(mut_off);             % Index of population/genotype that will produce mutants
                if  ~isempty(mi)
                    C        = C-mut_off;
                    C        = [C; ones(sum(mut_off),1)]; % Founding of mutant population
                    xe_m     = [];
                    xh_m     = [];
                    for pi   = mi'                        % Mutation process
                        meh  = mvnrnd([xe(pi),xh(pi)],sigma,mut_off(pi));
                        meh  = min(max(meh,0),1);
                        xe_m = [xe_m; meh(:,1)];
                        xh_m = [xh_m; meh(:,2)];
                    end
                    update_traits
                end
                
                
                %%% EXTINCTION %%%
                
                if  any(C==0)
                    ei     = find(~C);
                    C(ei)  = [];
                    xe(ei) = [];
                    xh(ei) = [];
                    e1(ei) = [];
                    e2(ei) = [];
                    h1(ei) = [];
                    h2(ei) = [];
                end
                
            end
            
            %%% BRANCHING TEST %%%
            spec1 = xe <= 0.4 & xh <=0.4;
            spec2 = xe >= 0.6 & xh >= 0.6;
            
            popspec1  = sum(C.*(spec1));
            popspec2  = sum(C.*(spec2));
            popTotal  = sum(C);
            
            
            if popspec1/popTotal > 0.35 & popspec2/popTotal > 0.35 % Branching
                V = [V;1];
            elseif  popspec1/popTotal > 0.90                       % Repelling
                V = [V;-1];
            elseif  isempty(C)                                     % Non-viable
                V = [V;0];
            else                                                   % CSS
                V = [V;-2];
            end
        end
        
        Data = [Data V]
        save('Data.mat', 'Data');
        
    end
    
end


%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%
    function update_traits
        xe = [xe; xe_m]; % Updating population's gene pool in locus e
        xh = [xh; xh_m]; % Updating population's gene pool in locus h
        e1 = [e1; 1 *     ((1 - exp(1) .^ (-te .* (1 - xe_m))) ./ (1 - exp(1) .^ (-te)))]; % Updating trait e1
        e2 = [e2; 1 *     ((1 - exp(1) .^ (-te .*      xe_m))  ./ (1 - exp(1) .^ (-te)))]; % Updating trait e2
        h1 = [h1; 1 - 1 * ((1 - exp(1) .^ (-th .* (1 - xh_m))) ./ (1 - exp(1) .^ (-th)))]; % Updating trait h1
        h2 = [h2; 1 - 1 * ((1 - exp(1) .^ (-th .*      xh_m))  ./ (1 - exp(1) .^ (-th)))]; % Updating trait h2
    end

fprintf('Simulation complete')
end