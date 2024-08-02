January 2023

## Black and Scholes equation
This project is a C++ implementation of finite difference methods to solve the parabolic partial differential equation (PDE) of the Black-Scholes model for pricing European call options.
The code includes functionality for option pricing, handling boundary conditions, and performing error analysis compared to analytical solutions.

### Project Structure
The project is organized as follows:
- **BS.h**: header file defining the BS class and its methods;
- **BS.cpp**: implementation of BS class methods;
- **main.cpp**: contains the main function and example usage.

### Boundary Conditions
A notable feature of this project is the handling of boundary conditions, which are not fixed and can move through a linear interpolation with a fixed point approximating infinity.

### Report
For an in-depth look at the results and methodologies used in this project, refer to the presentation **Black-Scholes.pdf** (in Italian).
