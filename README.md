# C-RANSAC
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Centralized RANSAC Based 3D Registration with Fast Convergence and High Accuracy (2023) by Kuo-Liang Chung and Wei-Tai Chang.  

<div align=center>
<img src="https://github.com/ivpml84079/C-RANSAC/blob/main/Fig/Example.png">
</div>

The left figure is the testing point cloud "Redkitchen" selected from the [3DMATCH](https://3dmatch.cs.princeton.edu/) dataset, the green lines are inliers between the source point cloud and the target point cloud, and the red lines are outliers. The right one is the registration result using our method.

## Acknowledgments
The programing implementation of our method is based on:

[H. Yang](http://hankyang.mit.edu/), [J. Shi](http://jingnanshi.com/), and [L. Carlone](http://lucacarlone.mit.edu/), “TEASER: Fast and Certifiable Point Cloud Registration,” *IEEE Trans. Robotics*, Dec. 2020. ([IEEE](https://ieeexplore.ieee.org/document/9286491)) ([github](https://github.com/MIT-SPARK/TEASER-plusplus))

We appreciate the authors for sharing their codes.

## TEASER++ Modified

We have made modifications to the TEASER++ library to showcase our method. The modified code is located in the "registration.cc" file. Additionally, we utilize and modify the "teaser_cpp_ply.cc" file for conducting experiments.

## Usage
Compiling the codes of our method in the project directory
```
sudo apt install cmake libeigen3-dev libboost-all-dev
git clone https://github.com/ivpml84079/C-RANSAC.git
cd C-RANSAC/TEASER-plusplus && mkdir build && cd build
sudo cmake .. && sudo make && sudo make install && sudo ldconfig
cd .. && cd examples/teaser_cpp_ply && mkdir build && cd build
sudo cmake .. && sudo make
```
using the execution code to estimate the registration solution
```
./teaser_cpp_ply {PATH_TO_POINT_CLOUD}
```
## Testing enviroment
* Windows Subsystem Linux (Ubuntu 20.04.6 LTS)
* ISO C++ 14

## Contact
If you have any questions, please email us via   
Wei-Tai Chang: m11115007@mail.ntust.edu.tw  
Kuo-Liang Chung: klchung01@gmail.com