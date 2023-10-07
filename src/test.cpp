#include <eigen3/Eigen/Dense>
#include <iostream>
int main()
{
    Eigen::Vector4d vec = {1,2,3,4};
    Eigen::seqN(0,3);
    std::cout << vec.norm() << '\n';
}