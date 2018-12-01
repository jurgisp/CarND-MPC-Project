# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

---

## Project Writeup

### MPC equations

I use 4 variables to model the sate:
* `x[t], y[t]` - coordinates (in car-relative frame)
* `psi[t]` - angle
* `v[t]` - velocity

And two actuator variables:
* `delta[t]` - turn angle
* `a[t]` - acceleration

The update equations are then:
```
x[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
y[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
psi[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
v[t] = v[t-1] + a[t-1] * dt
```

Diffrently from the examples in the class, I don't model CTE and Psi_error as explicit state variables, because there is no need - they don't have independent update equations, and are just intermediate variables when calculating the MPC cost.

### MPC cost function

The cost function is the following:
```
MPC_cost = 
  sum( (f(x[t]) - y[t])^2 ) +           // Cross-track error
  sum( (atan(f'(x[t])) - psi[t])^2 ) +  // Psi error
  sum( (v[t] - v_ref)^2 ) +             // Velocity error
  sum( (delta[t])^2 ) +                 // Minimize delta
  sum( (a[t])^2 ) +                     // Minimize acceleration
  sum( (delta[t] - delta[t-1])^2 ) +    // Minimize change in delta
  sum( (a[t] - a[t-1])^2 )              // Minimize change in acceleration
```

### MPC initial conditions

Initial conditions for state is trivial, because the car starts at (0, 0, 0) in car-relative coordinates:
* `x[0] = y[0] = psi[0] = 0`
* `v[0] = v_current`

Additionally we restrict the initial delta and a to be the last action:
* `delta[0] = delta_last`
* `a[0] = a_last`

Fixing action at `t=0` solves two problems simultaneously:

* It allows for the smooth action change term `sum((delta[t] - delta[t-1])^2)` to function properly. If we don't fix `delta[0]`, then we don't a way to ensure that the next `delta` action will be similar to the last.

* It *automatically* models the latency of one step, because the action that we return to the actuators are `(delta[1], a[1])`, which is the action that is appropriate `dt` time later.

### Handling latency

We solve the latency by fixing the action at the first step to be whatever it is. That way MPC itself assumes that we can not change the action for the first step (`dt` time), and only adjust the action at `(delta[1], a[1])`, which is `dt` time later.

Thus our choice of `dt` is correlated with the latency - we set it to be
* `dt = 150ms`

Which is, from measurement, the average time between actuator updates. That consists of `100ms` artificial lag, plus `~50ms` lag to solve MPC.

### Choice of (N, dt)

Our choice of `dt` is constrained by the latency modelling to be `150ms`, as described above.

That leaves us the choice of `N`, which is simply the tradeoff between time horizon (the bigger the better), and MPC solver latency (the bigger the worse). We choose `N=10` because that seems enough.

### Acceleration vs Throttle

There is a problem that MPC models actuator as acceleration `a[t]`, but the input that the car accepts is actually `throttle in [-1, 1]`. I found from experiments that:
* `throttle=1` results in `a=5m/s^2` at `v=0`
* For higher `v` same throttle results in smaller acceleration, up to a point where fixed `throttle` provides no acceleration, just keeps velocity constant.
* I didn't work out the exact functional correspondence `a(v,throttle)`, though it would be possible, and MPC could model that.

Without modelling the throttle properly, it is actually misleading to have `a[t]` as actuator variable, because we can't actually control the car correctly according to the MPC output.

In the end I worked around this by forcing MPC to assume no acceleration, and outputting constant `throttle=0.3` to the car. The car then accelerates to about `v~25mph`, at which point acceleration goes to zero, and the MPC model becomes correct, just by modelling `delta[t]`.

To improve this part, we would need to have `throttle[t]` in the MPC model as actuator variable, and derive the correct state update equations for that.

### Fitting waypoint polynomial

Polynomial `y=f[x]` is derived straighforwardly by:
1. Transforming global waypoint coordinates to car-local coordinates, using the global car reference coordinates.
1. Fitting the polynomial to car-relative waypoint coordinates

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.

* **Ipopt and CppAD:** Please refer to [this document](https://github.com/udacity/CarND-MPC-Project/blob/master/install_Ipopt_CppAD.md) for installation instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.

## Tips

1. It's recommended to test the MPC on basic examples to see if your implementation behaves as desired. One possible example
is the vehicle starting offset of a straight line (reference). If the MPC implementation is correct, after some number of timesteps
(not too many) it should find and track the reference line.
2. The `lake_track_waypoints.csv` file has the waypoints of the lake track. You could use this to fit polynomials and points and see of how well your model tracks curve. NOTE: This file might be not completely in sync with the simulator so your solution should NOT depend on it.
3. For visualization this C++ [matplotlib wrapper](https://github.com/lava/matplotlib-cpp) could be helpful.)
4.  Tips for setting up your environment are available [here](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/f758c44c-5e40-4e01-93b5-1a82aa4e044f/concepts/23d376c7-0195-4276-bdf0-e02f1f3c665d)
5. **VM Latency:** Some students have reported differences in behavior using VM's ostensibly a result of latency.  Please let us know if issues arise as a result of a VM environment.

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!

More information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/f1820894-8322-4bb3-81aa-b26b3c6dcbaf/lessons/b1ff3be0-c904-438e-aad3-2b5379f0e0c3/concepts/1a2255a0-e23c-44cf-8d41-39b8a3c8264a)
for instructions and the project rubric.

## Hints!

* You don't have to follow this directory structure, but if you do, your work
  will span all of the .cpp files here. Keep an eye out for TODOs.

## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to we ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./

## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).
