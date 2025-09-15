%%Monte Carlo Simulation of pi with Precision Control
%task: Estimate pi using Monte Carlo simulation until a specified precision
%(in significant figures) is reached. Repeat this for precision levels
%1 through 10 and analyze the growth of required iterations.
%
%Description:
%-For each precision level (significant figures), random points are
%sampled uniformly from the unit square [0,1] x [0,1].
%-The ratio of points falling inside a circle of radius 0.5 (centered
%at (0.5,0.5)) is used to estimate pi.
%-The simulation continues until the Monte Carlo estimate of pi matches
%a stable target value computed via the Gauss–Legendre algorithm,
%rounded to the requested number of significant figures.
%-This ensures we stop without using MATLAB's built-in pi constant.
%
%Inputs:
%None (all parameters are set internally in the script).
%
%Outputs:
%N - Vector of iteration counts required at each precision.
%estimatedPiVec - Vector of final Monte Carlo pi estimates.
%sigFigVec - Vector of significant figure levels used.
%
%Figures:
%A semilog plot (iterations vs. significant figures).
%
%Dependencies:
%roundSig(x,s) - Helper function to round to s significant figures.
%estimatePi(s) - Helper function using Gauss–Legendre to compute pi
%to s significant figures (for comparison target).
%
%Notes:
%Significant figures tested are 1 through 10.
%Using higher than 10 significant figures may be expensive computationally.
%
%References:
%Gauss–Legendre algorithm:
%https://en.wikipedia.org/wiki/Gauss–Legendre_algorithm

N = []; %initialize empty vector N to store how many points were used
estimatedPiVec = []; %initialize empty vector to keep the estimated pi per iteration
sigFigVec = []; %initialize empty vector to keep the significant figures used

%we take 1 through 10 significant figures. Getting 11 is too computationally
%costly
for sigFigs = 1:10
    %initialize for the while loop
    stable = false;
    %uses helper function to get a pi estimate, not real pi
    %we make sure this is done before the while loop, not inside of it!!!
    targetRounded = estimatePi(sigFigs);
    %initialize for counting
    totalPoints = 0;
    insidePoints = 0;

    while ~stable
        %similar to task 1
        totalPoints = totalPoints + 1;
        x = rand(1);
        y = rand(1);
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2);
    
        if dist <= 0.5
            insidePoints = insidePoints + 1;
        end
    
        piEstimate = 4 * (insidePoints / totalPoints);

        %check if we're at our desired accuracy
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

%another helper function to estimate pi
function piEstimate = estimatePi(sigFig)
    %estimate pi to a certain number of significant figures using
    %gauss-legendre
    %see https://en.wikipedia.org/wiki/Gauss–Legendre_algorithm for more
    %information

    a = 1.0;
    b = 1 / sqrt(2);
    t = 0.25;
    p2 = 1.0;

    prevRounded = NaN; %initialize prevRounded
    stable = false; %for the while loop

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
%use log scaling since we explode fast
semilogy(sigFigVec, N, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Significant Figures');
ylabel('Iterations Needed (log scale)');
title('Monte Carlo Estimiating pi: Iterations vs Significant Figures');
grid on;