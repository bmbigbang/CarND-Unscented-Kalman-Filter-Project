# Unscented Kalman Filter Project

The program is capable of combining multidimensional sensor data
together using first order Taylor expansion term of the equation models
of motion for each dimension, while scaling the influence of each sensor
 by its uncertainty. This method is based on a Kalman Filter
  [E. A. Wan et al.](https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf)
  and is more efficient than direct derivative space conversion of the
  equation models, especially in higher dimensional sensors.
  
Once run, the data is read from sample exports of LIDAR and RADAR sensors
and given the initial sensor error values, a map of locations is generated
where at each point the uncertainty of either sensor could vary depending
on environmental causes. For example, the LIDAR might have reduced vision in a construction site due to physical obstacles, and return a high uncertainty value or deviate substantially from the current result. In this case the 
Filter will automatically take more influence from the RADAR measurement
and return a value closer to the current belief as a result.

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`


