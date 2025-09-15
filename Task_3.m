function piEstimate = Task_3(sigFigs)
    
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


%-----all copied from Task_2.m-----
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