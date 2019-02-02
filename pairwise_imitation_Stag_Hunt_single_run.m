clear all

%%%%%%%%%%%
%% SETUP %%
%%%%%%%%%%%

N=20; % Population size 
p_revision=1/2; % Probability of action revision
p_error=0.05;  % Probability of mutation
N_rounds=50000; % Number of rounds
k=N-1; % Required for looping over P


% Players are repeatedly and randomly matched in pairs to play 
% a Stag Hunt stage game
% Stag=1, Hare=0

%  MOVES    PAYOFFS
% P1   P2   P1   P2
% A    B    C    D
% 1    1    10   10
% 1    0    0    3
% 0    1    3    0
% 0    0    3    3


% The evolution of play under Revision protocol 1 (PRI) is described by the
% array Action1

Action1=round(rand(1,N)); % Initial state
% Random initial condition: Action1=round(rand(1,N));
% Everyone plays Hare at round 1: Action1 = zeros(1,N);
% Everyone plays Stag at round 1: Action1 = ones(1,N);

round=1; % set first round of play

% The agents' payoffs under Revision protocol 1 (PRI) are described by the
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
% (PII) are described by the arrays Action2 and Payoff2, respectively
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
       
    % Every player has a probability p_revision of being selected for
    % action revision. To impose the revising agents to be the same under
    % the two protocols being considered, set revision2=revision1
    revision1=rand(1,N)<p_revision;
    revision2=rand(1,N)<p_revision;
    %revision2=revision1;
    
    % In every round, Player i plays Action(round,i) and receives a
    % payoff of Payoff(round,i).

    
    %%%%%%%% REVISION PROTOCOL 1: PAIRWISE RANDOM IMITATION (PRI) %%%%%%%%
  
    for i=1:N
        % Check whether agent i has been given a revision opportunity
        if revision1(i)==true % If yes,...
        other_players=[1:i-1,i+1:N]; % ...consider the set P\{i}...
        reference=other_players(ceil(rand()*k)); % ...and draw 
            % one reference at random from P\{i}.
            % Check whether the payoff earned by i's reference in the
            % previous round is higher than the payoff earned by i.
            %If yes, then i copies their reference
            if Payoff1(round-1,reference)>Payoff1(round-1,i)
            Action1(round,i)=Action1(round-1,reference);
            end
            % If agent i makes a mistake, they selects the action 
            % that is NOT prescribed by the revision protocol
             if rand()<p_error
                Action1(round,i)=abs(Action1(round,i)-1);
                % Action1(round,i)= randi([0,1]); 
                % Select if mistakes consist in choosing one action 
                % at random
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

    
    %%%% REVISION PROTOCOL 2: PAIRWISE INTERACTION-AND-IMITATION (PII) %%%%

    for i=1:N
        % Check whether agent i has been given a revision opportunity
        if revision2(i)==true % If yes,...
        opponent=mod(i,2)*(2)+i-1; % ...consider agent i's last opponent
            % Check whether the payoff earned by i's opponent at time
            % round-1 is higher than the payoff earned by i. 
            % If yes, then i copies their opponent
            if Payoff2(round-1,opponent)>Payoff2(round-1,i)
            Action2(round,i)=Action2(round-1,opponent);
            end
            % When an agent makes a mistake, they selects the action
            % that is NOT prescribed by the revision protocol
            if rand()<p_error
                Action2(round,i)=abs(Action2(round,i)-1);
                % Action2(round,i)= randi([0,1]); 
                % Select if mistakes consist in choosing one action 
                % at random
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EVOLUTION OF COOPERATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If Player i played Stag at time t under protocolx, we have:
% Actionx(i,t)=1
% So the proportions of players that played Stag at round t are:
%p1_Stag(t)=sum(Action1(t,:))/N;
%p2_Stag(t)=sum(Action2(t,:))/N;

for t=1:N_rounds
    p1_Stag(t)=sum(Action1(t,:))/N;
    p2_Stag(t)=sum(Action2(t,:))/N;
end


%%%%%%%%%%
%% PLOT %%
%%%%%%%%%%

%figure('Name','% of agents playing Stag');
%hold on
%plot(p1_Stag);
%plot(p2_Stag);
%legend PRI PII
%hold off

figure()
set(gcf,'color','w')
plotPRI=subplot(2,1,1)
plot(p1_Stag)
title('{\normalfont(PRI)}','Interpreter','latex', 'fontsize', 12)  
xlabel('Time') 
ylabel({'Proportion of stag hunters'})
plotPII=subplot(2,1,2)
plot(p2_Stag,'color',[0.8500 0.3250 0.0980])
title('{\normalfont(PII)}','Interpreter','latex','fontsize', 12)  
xlabel('Time') 
ylabel({'Proportion of stag hunters'})
axis([plotPRI plotPII],[0 N_rounds 0 1])



