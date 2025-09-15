# AMS595-Project-1
Monte Carlo Estimation of pi

This project explores different ways to approximate pi using Monte Carlo simulation in MATLAB.  
The general idea: randomly sample points in the unit square [0,1] × [0,1], count how many fall inside a circle of radius 0.5 centered at (0.5,0.5), and use this ratio to approximate π.  

Three tasks were completed:

Task 1: Error and Runtime Analysis (`task_1.m`)
Goal: Estimate π for increasing sample sizes `N` and analyze accuracy and runtime.  
Method:  
  Runs the simulation for `N = 100` to `100,000`.  
  Computes estimated π, absolute error from MATLAB’s `pi`, and elapsed time.  

Outputs: 
  `N`, `estimatedPiVec`, `errorVec`, `timeElapsed` in the workspace.  
  Two plots:
  Error (blue) and estimated π (red) vs. number of points.  
  Error vs. runtime.  

-----

Task 2: Precision-Based Stopping Rule (`task_2.m`)
Goal: Run the Monte Carlo simulation until π is estimated to a given number of significant figures (1–10).  
Method:
  Uses the Gauss–Legendre algorithm (via helper function `estimatePi`) to get a stable target π value for comparison.  
  Stops simulation once the Monte Carlo estimate matches the target at the specified precision.  
Outputs:  
  Vectors `N`, `estimatedPiVec`, and `sigFigVec`.  
  A semilog plot of iterations vs. significant figures.  

-----

Task 3: Precision with Visualization (`Task_3.m`)
Goal: Stop the Monte Carlo simulation once π is accurate to a user-specified number of significant figures, **and visualize the points**.  
Method: 
  User specifies precision level (`sigFigs`) as input.  
  Random points are plotted:  
    Red = inside the circle  
    Blue = outside the circle  
  Simulation continues until estimate stabilizes to the desired precision.  
Outputs:
  Returns the final Monte Carlo π estimate as `piEstimate`.  
  Figure showing points, circle boundary, and a results box (π rounded, unrounded estimate, and number of points).  
  Results printed to the Command Window.  

Example usage for Task 3:
```matlab
value = Task_3(3);   % Estimate pi to 3 significant figures
```

References:
Gauss–Legendre algorithm - https://en.wikipedia.org/wiki/Gauss–Legendre_algorithm
