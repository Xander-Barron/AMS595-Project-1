%%Monte Carlo simulation but stopping at certain precision levels
%We will be using significant figures as the precision levels
%For example pi=3 is 1 significant figure
%pi=3.1415 is 5 significant figures

N = []; %initialize empty vector N to store how many points were used
estimatedPiVec = []; %initialize empty vector to keep the estimated pi per iteration
sigFigVec = [];

for sigFigs = 1:10

    stable = false;
    targetRounded = estimatePi(sigFigs);
    totalPoints = 0;
    insidePoints = 0;

    while ~stable
        totalPoints = totalPoints + 1;
        x = rand(1);
        y = rand(1);
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2);
    
        if dist <= 0.5
            insidePoints = insidePoints + 1;
        end
    
        piEstimate = 4 * (insidePoints / totalPoints);
        
        if roundSig(piEstimate, sigFigs) == targetRounded
            N(end+1) = totalPoints;
            estimatedPiVec(end+1) = piEstimate;
            sigFigVec(end+1) = sigFigs;
            stable = true;
        end
    end
end


%helper function to round to the nearest signifcant figure
function y = roundSig(x, s)
    if x == 0
        %special case, shouldn't be used but just in case
        y = 0;
        return
    end
    %find order of magnitude of x
    k = floor(log10(abs(x)));
    %scaling x, rounding it then scaling it back up
    y = round(x / 10^k, s-1) * 10^k;
end

function piEstimate = estimatePi(sigFig)
    %estimate pi to a certain number of significant figures using
    %gauss-legendre

    a = 1.0;
    b = 1 / sqrt(2);
    t = 0.25;
    p2 = 1.0;

    prevRounded = NaN; %initialize prevRounded
    stable = false;

    while ~(stable)
        a1 = (a + b) / 2;
        b1 = sqrt(a * b);
        t1 = t - p2 * (a - a1)^2;
        p2 = 2 * p2;
        pi_k = ((a1 + b1)^2) / (4 * t1);

        currRounded = roundSig(pi_k, sigFig);

        %check if we are stable or not
        if ~isnan(prevRounded) && currRounded == prevRounded
            piEstimate = currRounded;
            stable = true;
        else
            prevRounded = currRounded;
        end

        %for next iteration
        a = a1;
        b = b1;
        t = t1;
    end
end

figure;
semilogy(sigFigVec, N, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Significant Figures');
ylabel('Iterations Needed (log scale)');
title('Monte Carlo Estimiating pi: Iterations vs Significant Figures');
grid on;