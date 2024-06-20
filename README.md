# Numeric-Differential-Equation-Solvers
This is a compilation of many differential equation solvers I have been making in MATLAB. 
I'm going to explain the context behind these programs.
This was started after watching a video on finding optimal paths of a spacecraft using some very confusing and complex means. I figured, "I have no idea what this guy just explained but it was very cool. I want to do it, but better." And so on the road to doing it, but better I realized how tough the challenge was. I had to learn completely new ways of numerically solving DE's. And I'm almost there!

In order to solve the problem I realized I needed to look for an optimal "path". A.K.A a function of some sort. Minimizing functions? That sounds exactly like the calculus of variations! After some math work with the Euler-Lagrange Equation I realized that the partial differential equation was a complete nightmare. Especially if I wanted to simulate the entire solar system. Each celestial body is another nightmare term in the equation. Thus, I knew I had to solve it numerically.

With any complex problem, the best thing you can do is attempt a similar, yet much easier form of the problem: solving 1 variable ODE's with boundary conditions. Which led me to the finite difference method. I found it was probably the best solution for the problem I was tackling. The finite difference method creates a system of nonlinear equations of which I have to solve. Which brings me to the first program I created.

NewtonRaphson1SecondOrderVarDESolver uses Newton-Raphson's method of solving nonlinear systems of equations to solve for the Y values that most closely satisfy the systems of equations, and as such, the Differential Equation. The main problem I ran into was how to get all of the derivative terms for the jacobian. Symbolically getting them all was out of the question. I realized I could just numerically get each derivative by using the literal definition of the derivative. Thanks Calc 1!
Anyway, this method worked pretty well for 1 variable Differential Equations. But I realized that it just wasn't as fast as I wanted. Not to mention I needed to add another variable and two more orders to reach the fourth order PDE that I wanted to solve.

Thus, we are brought to Broyden's method. This method is essentially Newton-Raphson method but instead of doing an inverse matrix operation every iteration, you approximate the inverse matrix each iteration. This is much faster.
This brings us to Broydens1VarSecondOrderDESolver.

Finally, I made the two variable version of Broydens method which brings us to:
Broydens2VarSecondOrderDESolver. From here, I need to add the functionality to include the fourth time derivative of x and y and write it in the context of the initial problem I am trying to solve.

As of 6/20/2024, I have finished the 2 Variable, fourth order DE solver using Broyden's method.

Broydens2VarFourthOrderDESolver.m and GreaterEfficiencyBroydensFourthOrder.m are both finished versions of the 2 var, 4th order solvers. (GreaterEfficiencyBroydensFourthOrder.m is much better though). The main struggle I had finishing this project was the PDE I was using was extremely complex. It was to the point that no matter what parameters I used, it just never converged and just created awful, inefficient graphs. Fortunately, I made a great simplification which doesn't technically truly solve for the lowest acceleration over time (or Delta-V), it does meaningfully approximate and function very similar to the original equation. (Instead of sqrt(Px^2+Py^2), it's just (Px^2+Py^2)). Removing the square root made the solver much much better. It also now works with the MOON! The moon moves over time, which is the main benefit to this program. You can add as many moving celestial bodies as needed and it will take into account their motion as the vehicle itself moves. Keep in mind that this isn't a program that calculates the trajectory based on initial conditions, it uses variational calculus to determine the optimal path function. So it feels really awesome that it works.

