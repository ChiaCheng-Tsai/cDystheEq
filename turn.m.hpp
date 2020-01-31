// Automatically translated using m2cpp 2.0 on 2020-01-31 15:49:19

#ifndef TURN_M_HPP
#define TURN_M_HPP

#include <armadillo>
using namespace arma ;

cx_vec turn(cx_vec w, double dt, cx_vec K2i, uvec j)
{
  cx_vec v ;
  v = arma::exp(K2i*dt)%w ;
  v(j-1).fill(0) ;
  return v ;
}
#endif