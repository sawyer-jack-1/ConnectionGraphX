
set(0,'DefaultFigureWindowStyle','docked')
clear all;
close all;
clc;

%%

CYCLE_LENGTH = 10;
DIM_CONNECTION = 5;
NUMBER_OF_TEST_TRIALS = 1e3;

%%

%Build a trivial connection graph on a cycle
G = GraphX.cycleGraphX(CYCLE_LENGTH);
G = ConnectionGraphX(G,DIM_CONNECTION);

figure;
plot(G.graphObject);
title("Cycle graph G");

%Change a single signature to make it unbalanced

%G = G.setEdgeConnection(1,2,flip(eye(DIM_CONNECTION)));
%randSig = flip(eye(DIM_CONNECTION));

randSig = ConnectionGraphX.getRandomSOMatrix(DIM_CONNECTION);
G = G.setEdgeConnection(1,2, randSig);

%%

%Instantiate a matrix to store the connection resistances once we calculate
%them
trialResistances = zeros(CYCLE_LENGTH,CYCLE_LENGTH, NUMBER_OF_TEST_TRIALS);

f = waitbar(0, "Please wait...");

%The definition of connection resistance may depend on some x\in S^{d-1}
%So we sample a bunch of them uniformly and calculate the corresponding
%resistance, then plot a histogram.
for k=1:NUMBER_OF_TEST_TRIALS
    waitbar(k / NUMBER_OF_TEST_TRIALS, f, "Sampling connection resistance w.r.t. random vectors...");
    R = zeros(CYCLE_LENGTH, CYCLE_LENGTH);
    x = randn(DIM_CONNECTION,1);
    x = x / norm(x);

    for i = 1:CYCLE_LENGTH
        for j = 1:CYCLE_LENGTH
               
            if G.adjacencyMatrix(i,j) == 1
        
                %Build the m^x_{ij} vector, analogue of e_{i}-e_{j}
                d = DIM_CONNECTION;
                m = zeros(CYCLE_LENGTH * d, 1);
                m((d*(i-1) + 1):(d * i),:) = x;
                m((d*(j-1) + 1):(d * j),:) = - G.connectionAdjacency((d*(j-1) + 1):(d * j),(d*(i-1) + 1):(d * i)).' * x;

                %The action of L^+ * m is replicated by looking for least
                %square solution to Lx = m. Since (m perp ker(L)) this will
                %be equivalent to L^+ * m and faster computationally.
                R(i,j) = (m.' * lsqminnorm(G.connectionLaplacian, m));
                
            end
    
        end
    end

    trialResistances(:,:,k) = R;
end    

%For posterity calculate the regular effective resistance.

normalEffRes = zeros(CYCLE_LENGTH, CYCLE_LENGTH);

for i = 1:CYCLE_LENGTH
    for j = 1:CYCLE_LENGTH
               
        if G.adjacencyMatrix(i,j) == 1
        
            A = eye(CYCLE_LENGTH);
            m = A(:,i) - A(:,j);

            %The action of L^+ * m is replicated by looking for least
            %square solution to Lx = m. Since (m perp ker(L)) this will
            %be equivalent to L^+ * m and faster computationally.
            normalEffRes(i,j) = (m.' * lsqminnorm(G.graphLaplacian, m));
                
        end
    
    end
end

waitbar(.2,f,"Calculating statistics & plots...");

fprintf("==> Connection matrix for edge {1,2}:\n");
randSig
fprintf("\n");

waitbar(.4,f,"Calculating statistics & plots...");

fprintf("==> Variance of R^sigma(i,j) across all trials:\n");
var(trialResistances,0,3)
fprintf("\n");

waitbar(.6,f,"Calculating statistics & plots...");

fprintf("==> Mean of R^sigma(i,j) across all trials:\n");
mean(trialResistances,3)
fprintf("\n");

waitbar(.8,f,"Calculating statistics & plots...");

fprintf("==> Regular combinatorial Laplacian effective resistance:\n");
normalEffRes
fprintf("\n")

waitbar(.9,f,"Calculating statistics & plots...");

figure;
histogram(trialResistances(1,2,:), 20, 'Normalization','pdf');
title("Resistances between modified edge \{1,2\}");

figure;
histogram(trialResistances(3,4,:), 20, 'Normalization','pdf');
title("Resistances between other edges");

close(f)



