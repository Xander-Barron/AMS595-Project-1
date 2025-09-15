%%Task 1
%Estimating pi using Monte Carlo Algorithm
%This script estimates the value of pi using a Monte Carlo method.
%Random points are generated in the unit square [0,1] x [0,1], and the
%ratio of points inside the unit circle (radius 0.5, centered at (0.5,0.5))
%is used to approximate pi.
%
%Description:
%The simulation is repeated for N values ranging from 100 to 100,000.
%For each N:
%A timer (tic/toc) is used to measure runtime.
%Pi is estimated from the fraction of points inside the circle.
%The absolute error from MATLAB's built-in pi is computed.
%Results are stored in vectors for later analysis.
%
%Variables:
%N - Vector of sample sizes used (number of random points).
%estimatedPiVec - Vector of estimated pi values for each N.
%errorVec - Vector of absolute errors |estimated pi - pi|.
%timeElapsed - Vector of elapsed times (seconds) per iteration.
%
%Figures:
%Scatter plot of estimated pi (red) and error (blue) vs. N.
%Scatter plot of error vs. time elapsed.


N = []; %initialize empty vector N to store how many points were used
errorVec = []; %initialize empty vector to store the errors we attain
timeElapsed = []; %initialize empty vector to measure the time elapsed per iteration
estimatedPiVec = []; %initialize empty vector to keep the estimated pi per iteration

for i = linspace(100, 100000, 999)
    tic; %start timer
    %repeat the Monte Carlo simulation for each value of N
    insideCircle = 0; %reset count for each iteration
    for j = 1:i
        x = rand(1); %random real number between 0 and 1
        y = rand(1);
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2); %if this distance is more than 0.5, we are outside the circle
        if dist <= 0.5
            insideCircle = insideCircle + 1; %since we're inside the circle, we add 1 to inside Circle
        end
    end
    estimatedPi = 4 * (insideCircle / i); %estimated pi from Monte Carlo
    piError = abs(estimatedPi - pi); %using absolute error

    N(end+1) = i; %store the current number of points used
    errorVec(end+1) = piError; %store the current error
    timeElapsed(end+1) = toc; %store the elapsed time for the current iteration
    estimatedPiVec(end+1) = estimatedPi; %store estimated pi calculations
end

%plot of estimated pi and its absolute error from true pi
h1 = scatter(N, errorVec, "b.");
hold on;
xlabel("Number of points thrown (N)"); 
%no need for y label as we're overlaying two plots
title("Monte Carlo pi Error");
h2 = scatter(N, estimatedPiVec, "r."); %estimated pi scatter overlayed
grid on;
%add legend
legend([h1 h2], {"Absolute error of estimated pi", "Estimated pi"}, "Location", "best")

%new figure to plot execution time against error
figure;
scatter(timeElapsed, errorVec, "b.");
%add labels and title to the second scatter plot
xlabel("Time Elapsed (s)");
ylabel("Error from pi");
title("Monte Carlo pi Error vs Time Elapsed");
grid on;

%% Task 2:
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

%% Task 3
%Task_3 is a Monte Carlo simulation to estimate pi until a desired
%precision level.
%This performs a Monte Carlo simulation to approximate pi by randomly
%sampling points in the [0,1]x[0,1] unit square. The simulation stops when
%the estimated pi matches the true value of pi rounded to the specified
%number of significant figures. The true value of pi is coputed using the
%Gauss-Legendre algorithm.
%
%Input: sigFigs - Positive integer spcifying the number of significant
%figures required
%
%Output: piEstimate - Final unrounded Monte Carlo estimate of pi when we
%match the sigFigs required
%
%Side Effects:
%-A figure is produced which shows red points that lay inside a circle with
%radius=0.5 centered at (0.5, 0.5); blue points that lay outside the circle
%but are still within [0,1]x[0,1]; and the circle boundary.
%-A text box is next to graph that displays the rounded value of pi, the
%unrounded value of pi and the final number of random points used to
%achieve those numbers.
%-Results are also printed in the command window
%
%Dependencies:
%roundSig(x,s) - Helper function to round to s significant figures.
%estimatePi(s) - Helper function using Gauss–Legendre to compute pi
%to s significant figures (for comparison target).
%
%Example:
%To estimate pi to 5 significant figures
%value = Task_3(5)
%
%References
%Gauss–Legendre algorithm:
%https://en.wikipedia.org/wiki/Gauss–Legendre_algorithm

function piEstimate = Task_3(sigFigs)

    %ensure sigFigs is >=1 and a number
    validateattributes(sigFigs, {'numeric'}, {'scalar', 'integer', '>=', 1});
    
    %get estimated pi from below helper function
    targetRounded = estimatePi(sigFigs); 

    N = 0; %initialize N
    insidePoints = 0; %initialize insidePoints

    %initialize vectors to store points inside and outside
    xIn = [];
    yIn = [];
    xOut = [];
    yOut = [];
    %for the while loop
    stable = false;

    while ~stable
        %similar to task 1 and 2
        N = N +1;
        x = rand;
        y = rand;

        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2);
    
        if dist <= 0.5
            insidePoints = insidePoints + 1;
            xIn(end + 1) = x;
            yIn(end + 1) = y;
        else
            xOut(end + 1) = x;
            yOut(end + 1) = y;
        end

        piEstimate = 4 * (insidePoints / N);

        %stop when we match the target
        if roundSig(piEstimate, sigFigs) == targetRounded
            stable = true; %break the while loop
        end
    end

    %similar to task 2
    roundedPi = roundSig(piEstimate, sigFigs);

    %set up plot
    figure;
    hold on;
    axis equal;
    box on;
    %ensure we only see [0,1]x[0,1]
    xlim([0 1]);
    ylim([0 1]);
    xlabel("x");
    ylabel("y");
    title(sprintf("Monte Carlo pi Estimation with %d sig figs)", sigFigs));
    grid on;

    %create the cricle outline
    theta = linspace(0, 2 * pi, 500); %theta as in angle
    plot(0.5 + 0.5 * cos(theta), 0.5 + 0.5 * sin(theta), "k-", "LineWidth", 2);


    %scatter points with the colors
    scatter(xIn, yIn, 4, "filled", "r"); %r for red
    scatter(xOut, yOut, 4, "filled", "b"); %b for blue
    %legend for the plot
    legend({"Circle boundary", "Inside circle (red)", "Outside circle (blue)"}, "Location", "southoutside")

    %annotation for the computed values of pi and the total number of
    %points used
    anno = sprintf('pi = %.*g (to %d sig figs) \nUnrounded: %.10f\nN = %d', sigFigs, roundedPi, sigFigs, piEstimate, N);
    annotation("textbox", [0.75, 0.25, 0.2, 0.2], "String", anno, "FitBoxToText", "on", "Interpreter", "none", "FontWeight", "bold", "Margin", 6);

    %command window output
    fprintf("Final pi (rounded to %d sig figs): %.*g\n", sigFigs, sigFigs, roundedPi);
    fprintf("Unrounded estimate: %.10f with N = %d\n", piEstimate, N);
end

%user prompt for task 3
sigFigFor3 = input('Enter desired number of significant figures (>=1): ');

value = Task_3(sigFigFor3);


%the below is commented because it already exists
% 
% %-----all copied from Task_2.m-----
% %helper function to round to the nearest signifcant figure
% function y = roundSig(x, s)
%     if x == 0
%         %special case, shouldn't be used but just in case
%         y = 0;
%         return
%     end
%     %find order of magnitude of x
%     k = floor(log10(abs(x)));
%     %scaling x, rounding it then scaling it back up
%     y = round(x / 10^k, s-1) * 10^k;
% end
% 
% %another helper function to estimate pi
% function piEstimate = estimatePi(sigFig)
%     %estimate pi to a certain number of significant figures using
%     %gauss-legendre
%     %see https://en.wikipedia.org/wiki/Gauss–Legendre_algorithm for more
%     %information
% 
%     a = 1.0;
%     b = 1 / sqrt(2);
%     t = 0.25;
%     p2 = 1.0;
% 
%     prevRounded = NaN; %initialize prevRounded
%     stable = false; %for the while loop
% 
%     while ~(stable)
%         a1 = (a + b) / 2;
%         b1 = sqrt(a * b);
%         t1 = t - p2 * (a - a1)^2;
%         p2 = 2 * p2;
%         pi_k = ((a1 + b1)^2) / (4 * t1);
% 
%         currRounded = roundSig(pi_k, sigFig);
% 
%         %check if we are stable or not
%         if ~isnan(prevRounded) && currRounded == prevRounded
%             piEstimate = currRounded;
%             stable = true;
%         else
%             prevRounded = currRounded;
%         end
% 
%         %for next iteration
%         a = a1;
%         b = b1;
%         t = t1;
%     end
% end