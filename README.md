# Numeric-Differential-Equation-Solvers
This is a compilation of many differential equation solvers I have been making in MATLAB. 
I'm going to explain the context behind these programs.
This was started after watching a video on finding optimal paths of a spacecraft using some very confusing and complex means. I figured, "I have no idea what this guy just explained but it was very cool. I want to do it, but better." And so on the road to doing it, but better I realized how tough the challenge was. I had to learn completely new ways of numerically solving DE's. And I'm almost there!

In order to solve the problem I realized I needed to look for an optimal "path". A.K.A a function of some sort. Minimizing functions? That sounds exactly like the calculus of variations! After some math work with the Euler-Lagrange Equation I realized that the partial differential equation was a complete nightmare. Especially if I wanted to simulate the entire solar system. Each celestial body is another nightmare term in the equation. Thus, I knew I had to solve it numerically.

With any complex problem, the best thing you can do is attempt a similar, yet much easier form of the problem: solving 1 variable ODE's with boundary conditions. Which led me to the finite difference method. I found it was probably the best solution for the problem I was tackling. The finite difference method creates a system of nonlinear equations of which I have to solve. Which brings me to the first program I created.

NewtonRaphson1SecondOrderVarDESolver uses Newton-Raphson's method of solving nonlinear systems of equations to solve for the Y values that most closely satisfy the systems of equations, and as such, the Differential Equation. The main problem I ran into was how to get all of the derivative terms for the jacobian. Symbolically getting them all was out of the question. I realized I could just numerically get each derivative by using the literal definition of the derivative. Thanks Calc 1!
Anyway, this method worked pretty well for 1 variable Differential Equations. But I realized that it just wasn't as fast as I wanted. Not to mention I needed to add another variable and two more orders to reach the fourth order PDE that I wanted to solve.

