// Automatically translated using m2cpp 2.0 on 2020-01-31 15:49:19

#ifndef RHS_M_HPP
#define RHS_M_HPP

#include "turn.m.hpp"
#include <armadillo>
using namespace arma ;

cx_vec RHS(cx_vec w, double dt, double e, uvec j, double g, vec k, vec op, cx_vec K2i)
{
  cx_vec Hv, rh, rhs, v, v_hat, va2, vx ;
  v_hat = turn(w, -dt, K2i, j) ;
  v = arma::ifft<cx_vec>(v_hat) ;
  vx = arma::ifft<cx_vec>(cx_double(0, 1)*k%v_hat) ;
  va2 = arma::conj(v)%v ;
  Hv = arma::ifft<cx_vec>(op%arma::fft(va2)) ;
  rh = arma::fft(-va2%v*cx_double(0, 1)-va2%vx*8.0*e*g-v%Hv*cx_double(0, 4)*e*g) ;
  rhs = turn(rh, dt, K2i, j) ;
  return rhs ;
}
#endif