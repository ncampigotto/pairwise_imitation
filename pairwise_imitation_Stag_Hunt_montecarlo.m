clear all
clc
tic

% !!! Statistics and Machine Learning toolbox required !!! %
% This version: v.1.0 - No parallel computing allowed


%%%%%%%%%%%
%% SETUP %%
%%%%%%%%%%%

% Players are repeatedly and randomly matched in pairs to play a stag hunt
% stage game
% ACTIONS: Stag=1, Hare=0

% ACTIONS   PAYOFFS
% P1   P2   P1   P2
% 1    1    10   10
% 1    0    0    3
% 0    1    3    0
% 0    0    3    3

N=20; % Population size
p_revision=0.5; % Probability of action revision
p_error=0.05 ;  % Probability of mutation
N_rounds=50000;    % Number of rounds
k=N-1;  % Required for looping over P
N_repetitions=1000; % Number of Monte Carlo repetitions
                    % For the typical run, set N_repetitions=1

% The evolution of play under Revision protocol 1 (PRI) is stored in the
% array Action1

% Start Monte Carlo experiment
for repetition=1:N_repetitions

% Initial state
Action1 = [zeros(1,N/2) ones(1,N/2)];
Action1 = Action1(randperm(numel(Action1))); 

% - Random initial condition: Action1=round(rand(1,N));
% - Everyone plays Hare at round 1: Action1 = zeros(1,N);
% - Everyone plays Stag at round 1: Action1 = ones(1,N);
% - 1/2 of the population plays Stag and 1/2 of the population plays Hare
%   at round 1: Action1 = [zeros(1,N/2) ones(1,N/2)];
%               Action1 = Action1(randperm(numel(Action1)));

round=1; % Set first round of play

% The agents' payoffs under Revision protocol 1 (PRI) are stored in the
% array Payoff1
% Remark: agent i (odd) plays against agent i+1 (even)
for i=1:2:N-1
    if(Action1(round,i)==1 && Action1(round,i+1)==1)
    Payoff1(round,i)=10;
    Payoff1(round,i+1)=10;
    end
    if(Action1(round,i)==0 && Action1(round,i+1)==0)
    Payoff1(round,i)=3;
    Payoff1(round,i+1)=3;
    end
    if(Action1(round,i)==0 && Action1(round,i+1)==1)
    Payoff1(round,i)=3;
    Payoff1(round,i+1)=0;
    end
    if(Action1(round,i)==1 && Action1(round,i+1)==0)
    Payoff1(round,i)=0;
    Payoff1(round,i+1)=3;
    end
end

% The evolution of play and the agents' payoffs under Revision protocol 2
% (PII) are stored in the arrays Action2 and Payoff2, respectively
Action2=Action1; 
Payoff2=Payoff1;


%%%%%%%%%%
%% LOOP %%
%%%%%%%%%%

for round=2:N_rounds
    
    % Initialise play:
    Action1(round,:)=Action1(round-1,:);
    Action2(round,:)=Action2(round-1,:);
    Payoff1(round,:)=Payoff1(round-1,:);
    Payoff2(round,:)=Payoff2(round-1,:);
       
    % Every player has a probability p_revision of receiving a revision
    % opportunity. To impose the revising agents to be the same under
    % the two protocols being considered, set revision2=revision1
    
    revision1=rand(1,N)<p_revision;
    revision2=rand(1,N)<p_revision;
    %revision2=revision1;
    
    % In every round, Player i plays Action(round,i) and receives a
    % payoff of Payoff(round,i).
        
    
    %%% REVISION PROTOCOL 1: PAIRWISE RANDOM IMITATION (PRI) %%%%
  
    for i=1:N
        % Check whether agent i has been given a revision opportunity. If
        % yes...
        if revision1(i)==true
        % ...consider the set P\{i}...
        other_players=[1:i-1,i+1:N];
        % ...and draw one reference at random from P\{i}
        reference=other_players(ceil(rand()*k));  
            % Check whether the payoff earned by i's reference in the
            % previous round is higher than the payoff earned by i. If yes,
            % then i copies her reference
            if Payoff1(round-1,reference)>Payoff1(round-1,i)
            Action1(round,i)=Action1(round-1,reference);
            end
            % If agent i makes a mistake, she selects the action that is
            % NOT prescribed by the revision protocol
             if rand()<p_error
                Action1(round,i)=abs(Action1(round,i)-1);
                % If, instead, mutants are assumed to choose one action at
                % random, set Action1(round,i)= randi([0,1]); 
             end
        end
    end

    % Randomly shuffle Action1/Payoff1 (new pairs are formed at random)
    lastPlayPayoff1=[Action1(end,:);Payoff1(end,:)];
    shuffled_lastPlayPayoff1 = lastPlayPayoff1(:,randperm(size(lastPlayPayoff1,2)));
    Action1(end,:)=shuffled_lastPlayPayoff1(end-1,:);
    Payoff1(end,:)=shuffled_lastPlayPayoff1(end,:);
    clear lastPlayPayoff1 shuffled_lastPlayPayoff1;
        
    % Compute payoffs
    for i=1:2:N-1
        if(Action1(round,i)==1 && Action1(round,i+1)==1)
        Payoff1(round,i)=10;
        Payoff1(round,i+1)=10;
        end
        if(Action1(round,i)==0 && Action1(round,i+1)==0)
        Payoff1(round,i)=3;
        Payoff1(round,i+1)=3;
        end
        if(Action1(round,i)==0 && Action1(round,i+1)==1)
        Payoff1(round,i)=3;
        Payoff1(round,i+1)=0;
        end
        if(Action1(round,i)==1 && Action1(round,i+1)==0)
        Payoff1(round,i)=0;
        Payoff1(round,i+1)=3;
        end
    end

    
    %%% REVISION PROTOCOL 2: PAIRWISE INTERACTION-AND-IMITATION (PII) %%%

    for i=1:N
        % Check whether agent i has been given a revision opportunity. If
        % yes...
        if revision2(i)==true
        % ...consider agent i's last opponent
        opponent=mod(i,2)*(2)+i-1;
            % Check whether the payoff earned by i's opponent at time
            % round-1 is higher than the payoff earned by i. If yes, then i
            % copies her opponent
            if Payoff2(round-1,opponent)>Payoff2(round-1,i)
            Action2(round,i)=Action2(round-1,opponent);
            end
            % If agent i makes a mistake, she selects the action that is
            % NOT prescribed by the revision protocol
            if rand()<p_error
                Action2(round,i)=abs(Action2(round,i)-1);
                % If, instead, mutants are assumed to choose one action at
                % random, set Action2(round,i)= randi([0,1]); 
            end
        end
    end
    
    % Randomly shuffle Action2/Payoff2 (new pairs are formed at random)
    lastPlayPayoff2=[Action2(end,:);Payoff2(end,:)];
    shuffled_lastPlayPayoff2 = lastPlayPayoff2(:,randperm(size(lastPlayPayoff2,2)));
    Action2(end,:)=shuffled_lastPlayPayoff2(end-1,:);
    Payoff2(end,:)=shuffled_lastPlayPayoff2(end,:);
    clear lastPlayPayoff2 shuffled_lastPlayPayoff2;
    
    % Compute payoffs
    for i=1:2:N-1
        if(Action2(round,i)==1 && Action2(round,i+1)==1)
        Payoff2(round,i)=10;
        Payoff2(round,i+1)=10;
        end
        if(Action2(round,i)==0 && Action2(round,i+1)==0)
        Payoff2(round,i)=3;
        Payoff2(round,i+1)=3;
        end
        if(Action2(round,i)==0 && Action2(round,i+1)==1)
        Payoff2(round,i)=3;
        Payoff2(round,i+1)=0;
        end
        if(Action2(round,i)==1 && Action2(round,i+1)==0)
        Payoff2(round,i)=0;
        Payoff2(round,i+1)=3;
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%
%% HISTORIES OF PLAY %%
%%%%%%%%%%%%%%%%%%%%%%%

% If Player i played Stag at time t under protocol a, we have:
% Actiona(i,t)=1.
% So the proportions of players that played Stag at round t are:
% p1_Stag_round(t)=sum(Action1(t,:))/N (under PRI)
% p2_Stag_round(t)=sum(Action2(t,:))/N (under PRI)

for t=1:N_rounds
    p1_Stag_round(t)=sum(Action1(t,:))/N;
    p2_Stag_round(t)=sum(Action2(t,:))/N;
end

% If Num_repetitions=1 and we want the typical history of play, enable the
% following. 

%figure()
%set(gcf,'color','w')
%plotPRI=subplot(2,1,1)
%plot(p1_Stag_round)
%title('{\normalfont(PRI)}','Interpreter','latex', 'fontsize', 12)  
%xlabel('Time') 
%ylabel({'Proportion of stag hunters'})
%plotPII=subplot(2,1,2)
%plot(p2_Stag_round,'color',[0.8500 0.3250 0.0980])
%title('{\normalfont(PII)}','Interpreter','latex','fontsize', 12)  
%xlabel('Time') 
%ylabel({'Proportion of stag hunters'})
%axis([plotPRI plotPII],[0 N_rounds 0 1])

% How does the proportion of stag hunters vary over time?

if repetition==1 % In the first repetition...
    p1_Stag=p1_Stag_round; % ...create new history of play
    p2_Stag=p2_Stag_round;
else % From repetition 2 onwards... 
    p1_Stag=[p1_Stag ; p1_Stag_round]; % ...attach the history from the last repetition
    p2_Stag=[p2_Stag ; p2_Stag_round]; % to the matrix of all histories
end

clear round % Avoids errors when defining Action1=round(rand(1,N))

end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AVERAGE PROPORTION OF STAG HUNTERS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute average proportion of stag hunters under PRI and PII
avg_p1_Stag=mean(p1_Stag);
avg_p2_Stag=mean(p2_Stag);

% Compute bootstrapped confidence intervals at the 95% level
CI_p1_Stag = bootci(10000, @mean, p1_Stag);
CI_p2_Stag = bootci(10000, @mean, p2_Stag);


%%%%%%%%%%%
%% PLOTS %%
%%%%%%%%%%%

figure()
set(gcf,'color','w')
plotPRI=subplot(2,1,1)
plot(CI_p1_Stag(1,:), '--', 'color', [.75 .75 .75])
hold on
plot(CI_p1_Stag(2,:), '--', 'color', [.75 .75 .75])
plot(avg_p1_Stag,'color',[0 0.4470 0.7410])
hold off
title('{\normalfont(PRI)}','Interpreter','latex', 'fontsize', 12)  
xlabel('Time') 
ylabel({'% of stag hunters'})
plotPII=subplot(2,1,2)
plot(CI_p2_Stag(1,:), '--', 'color', [.75 .75 .75])
hold on
plot(CI_p2_Stag(2,:), '--', 'color', [.75 .75 .75])
plot(avg_p2_Stag,'color',[0.8500 0.3250 0.0980])
hold off
title('{\normalfont(PII)}','Interpreter','latex','fontsize', 12)  
xlabel('Time') 
ylabel({'% of stag hunters'})
axis([plotPRI plotPII],[0 N_rounds 0 1])



