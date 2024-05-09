# The Poisson equation

 Here we consider the Poisson problem with some OpenMP parallelism:

$$
 \frac{{\partial^2 u}}{{\partial x^2}} + \frac{{\partial^2 u}}{{\partial y^2}} + \frac{{\partial^2 u}}{{\partial z^2}} = -f(x, y, z), \quad (x, y, z) \in \Omega
$$

where $u = u(x, y, z)$ is the solution, $f = f(x, y, z)$ is a source term, and $\Omega$ is the problem domain.

This code calculates the heat distribution in a cubic room with a radiator placed near a cold wall. Convection and similar effects are ignored. 

The domain $\Omega$ is defined

$$
\Omega = {(x, y, z) : |x| \leq 1, |y| \leq 1, |z| \leq 1} 
$$

We use Dirichlet boundary conditions defined below. 

$$
u(x, 1, z) = 20, u(x, -1, z) = 0 \qquad (|x| \leq 1), (|z| \leq 1)
$$

$$
u(1, y, z) = u(-1, y, z) = 20) \qquad (|y| \leq 1), (|z| \leq 1)
$$

$$
u(x, y, -1) = u(x, y, 1) = 20) \qquad (|x| \leq 1), (|y| \leq 1)
$$

Finally, the radiator is represented by the function:

$$
f(x, y, z) = \begin{cases} 200, & -1 \leq x \leq -\frac{3}{8}, \quad  -1 \leq y \leq -\frac{1}{2}, \quad -\frac{2}{3} \leq z \leq 0 \ \\
0, & \text{elsewhere} \end{cases}
$$

<img
  src="radiator.png"
  style="display: block; margin: 0 auto; max-width: 300px">

The problem is discretised on a  cubic $N \times N \times N$ grid, where $(u_{i,j,k})$ gives the solution at the point $(i,j,k)$. 

### Parallelism

OpenMP is implemented to include some lightweight parallelism. The code also produces a file `performance.csv` giving basic speedup information. 


***
### Credit

This solver was written by Elloïse Fangel-Lloyd, Dana Lüdemann, and Clemens Zengler.
