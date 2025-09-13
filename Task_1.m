%% Estimating pi using Monte Carlo algorithm

N = []; %initialize empty vector N to store how many points were used
errorVec = []; %initialize empty vector to store the errors we attain
timeElapsed = []; %initialize empty vector to measure the time elapsed per iteration
estimatedPiVec = []; %initialize empty vector to keep the estimated pi per iteration

for i = linspace(100, 100000, 999)

    tic; %start timer
    %repeat the Monte Carlo simulation for each value of N
    insideCircle = 0; %reset count for each iteration
    for j = 1:i
        x = rand(1);
        y = rand(1);
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2);
        if dist <= 0.5
            insideCircle = insideCircle + 1;
        end
    end
    estimatedPi = 4 * (insideCircle / i);
    piError = abs(estimatedPi - pi);

    N(end+1) = i; %store the current number of points used
    errorVec(end+1) = piError; %store the current error
    timeElapsed(end+1) = toc; %store the elapsed time for the current iteration
    estimatedPiVec(end+1) = estimatedPi;
end

scatter(N, errorVec, "b.");
hold on;
xlabel("Number of points thrown (N)");
ylabel("Error from pi");
title("Monte Carlo pi Error");
scatter(N, estimatedPiVec, "r.");
grid on;

figure;
scatter(timeElapsed, errorVec, "b.");
%add labels and title to the second scatter plot
xlabel("Time Elapsed (s)");
ylabel("Error from pi");
title("Monte Carlo pi Error vs Time Elapsed");
grid on;