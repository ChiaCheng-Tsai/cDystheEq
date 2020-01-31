// Automatically translated using m2cpp 2.0 on 2020-01-31 15:49:19

/*
%{
    A pseudo-spectral solver for Dysthe-Lo-Mei equation.
    See Equations (2.6) - (2.9) from Lo & Mei (1985, JFM)
    paper for more details.

    Here we test the code by simulating the NLS ground state evolution
    and deformation under higher order Dysthe terms. This solution is
    exact if the parameter e = 0.

    Copyright (C) 2020 Denys DUTYKH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

%}

%%% -------------------------------------------------- %%%
%%% The main file for pseudo-spectral Dysthe solver    %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, USMB           %%%
%%% E-mail: Denys.Dutykh@univ-smb.fr                   %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%
%%% Distributed under GNU General Public License       %%%
%%% -------------------------------------------------- %%%
*/
#include "turn.m.hpp"
#include "RHS.m.hpp"
#include "matlab2cpp.h"
#include <armadillo>
using namespace arma ;

int main(int argc, char** argv)
{
  cx_vec K2i, k1, k10, k11, k12, k13, k14, k15, k16, k2, k3, k4, k5, k6, k7, k8, k9, u, u_hat, v, w ;
  double C36, Tf, a, a10, a11, a12, a13, a15, a2, a3, a4, a5, a6, a7, a8, a9, b10_1, b10_6, b10_7, b10_8, b10_9, b11_1, b11_10, b11_7, b11_8, b11_9, b12_1, b12_10, b12_11, b12_6, b12_7, b12_8, b12_9, b13_1, b13_10, b13_11, b13_12, b13_6, b13_7, b13_8, b13_9, b14_1, b14_10, b14_11, b14_12, b14_13, b14_6, b14_7, b14_8, b14_9, b15_1, b15_10, b15_11, b15_12, b15_13, b15_6, b15_7, b15_8, b15_9, b16_1, b16_10, b16_11, b16_12, b16_13, b16_15, b16_6, b16_7, b16_8, b16_9, b21, b31, b32, b41, b43, b51, b53, b54, b61, b64, b65, b71, b74, b75, b76, b81, b86, b87, b91, b96, b97, b98, c1, c10, c11, c12, c13, c14, c8, c9, dt, dt2, dx, e, e1, e10, e11, e12, e13, e14, e15, e16, e8, e9, er1, er2, g, h, l, om, rho1, rho2, s6, t, tdraw, tloc, tol ;
  int N ;
  uvec j ;
  vec ONES, k, op, u0, x ;

  e = 0.05 ;
  g = 1.0 ;
  h = 10.0 ;
  N = 2048 ;
  l = 40.0 ;
  dx = 2*l/N ;
  x = arma::trans((m2cpp::fspan(1-N/2.0, 1, N/2.0)))*dx ;
  k = arma::trans(arma::join_rows(arma::join_rows(m2cpp::fspan(0, 1, N/2.0-1), m2cpp::srow<double>(0)), m2cpp::fspan(1-N/2.0, 1, -1)))*datum::pi/l ;
  op = -0.5*k%tanh(2.0*k*h) ;
  j = m2cpp::span<uvec>(N/4+2, 1, N/4*3) ;
  k(j-1).fill(0) ;
  K2i = arma::square(k)*(-cx_double(0, 1)*(cx_double) g*(cx_double) g) ;
  s6 = sqrt(6.0) ;
  c1 = 103.0/1680.0 ;
  c8 = -27.0/140.0 ;
  c9 = 76.0/105.0 ;
  c10 = -201.0/280.0 ;
  c11 = 1024.0/1365.0 ;
  c12 = 3.0/7280.0 ;
  c13 = 12.0/35.0 ;
  c14 = 9.0/280.0 ;
  a2 = 1.0/12.0 ;
  a3 = 1.0/9.0 ;
  a4 = 1.0/6.0 ;
  a5 = 2.0*(1.0+s6)/15.0 ;
  a6 = (6.0+s6)/15.0 ;
  a7 = (6.0-s6)/15.0 ;
  a8 = 2.0/3.0 ;
  a9 = 1.0/2.0 ;
  a10 = 1.0/3.0 ;
  a11 = 1.0/4.0 ;
  a12 = 4.0/3.0 ;
  a13 = 5.0/6.0 ;
  a15 = 1.0/6.0 ;
  b21 = 1.0/12.0 ;
  b31 = 1.0/27.0 ;
  b32 = 2.0/27.0 ;
  b41 = 1.0/24.0 ;
  b43 = 1.0/8.0 ;
  b51 = (4.0+94*s6)/375.0 ;
  b53 = -(94.0+84*s6)/125.0 ;
  b54 = (328.0+208*s6)/375.0 ;
  b61 = (9.0-s6)/150.0 ;
  b64 = (312.0+32*s6)/1425.0 ;
  b65 = (69.0+29*s6)/570.0 ;
  b71 = (927.0-347*s6)/1250.0 ;
  b74 = (-16248.0+7328*s6)/9375.0 ;
  b75 = (-489.0+179*s6)/3750.0 ;
  b76 = (14268.0-5798*s6)/9375.0 ;
  b81 = 2.0/27.0 ;
  b86 = (16.0-s6)/54.0 ;
  b87 = (16.0+s6)/54.0 ;
  b91 = 19.0/256.0 ;
  b96 = (118.0-23*s6)/512.0 ;
  b97 = (118.0+23*s6)/512.0 ;
  b98 = -9.0/256.0 ;
  b10_1 = 11.0/144.0 ;
  b10_6 = (266.0-s6)/864.0 ;
  b10_7 = (266.0+s6)/864.0 ;
  b10_8 = -1.0/16.0 ;
  b10_9 = -8.0/27.0 ;
  b11_1 = (5034.0-271*s6)/61440.0 ;
  b11_7 = (7859.0-1626*s6)/10240.0 ;
  b11_8 = (-2232.0+813*s6)/20480.0 ;
  b11_9 = (-594.0+271*s6)/960.0 ;
  b11_10 = (657.0-813*s6)/5120.0 ;
  b12_1 = (5996.0-3794*s6)/405.0 ;
  b12_6 = -(4342.0+338*s6)/9.0 ;
  b12_7 = (154922.0-40458*s6)/135.0 ;
  b12_8 = (-4176.0+3794*s6)/45.0 ;
  b12_9 = (-340864.0+242816*s6)/405.0 ;
  b12_10 = (26304.0-15176*s6)/45.0 ;
  b12_11 = -26624.0/81.0 ;
  b13_1 = (3793.0+2168*s6)/103680.0 ;
  b13_6 = (4042.0+2263*s6)/13824.0 ;
  b13_7 = (-231278.0+40717*s6)/69120.0 ;
  b13_8 = (7947.0-2168*s6)/11520.0 ;
  b13_9 = (1048.0-542*s6)/405.0 ;
  b13_10 = (-1383.0+542*s6)/720.0 ;
  b13_11 = 2624.0/1053.0 ;
  b13_12 = 3.0/1664.0 ;
  b14_1 = -137.0/1296.0 ;
  b14_6 = (5642.0-337*s6)/864.0 ;
  b14_7 = (5642.0+337*s6)/864.0 ;
  b14_8 = -299.0/48.0 ;
  b14_9 = 184.0/81.0 ;
  b14_10 = -44.0/9.0 ;
  b14_11 = -5120.0/1053.0 ;
  b14_12 = -11.0/468.0 ;
  b14_13 = 16.0/9.0 ;
  b15_1 = (33617.0-2168*s6)/518400.0 ;
  b15_6 = (-3846.0+31*s6)/13824.0 ;
  b15_7 = (155338.0-52807*s6)/345600.0 ;
  b15_8 = (-12537.0+2168*s6)/57600.0 ;
  b15_9 = (92.0+542*s6)/2025.0 ;
  b15_10 = -(1797.0+542*s6)/3600.0 ;
  b15_11 = 320.0/567.0 ;
  b15_12 = -1.0/1920.0 ;
  b15_13 = 4.0/105.0 ;
  b16_1 = -(36487.0+30352*s6)/279600.0 ;
  b16_6 = -(29666.0+4499*s6)/7456.0 ;
  b16_7 = (2779182.0-615973*s6)/186400.0 ;
  b16_8 = (-94329.0+91056*s6)/93200.0 ;
  b16_9 = (-232192.0+121408*s6)/17475.0 ;
  b16_10 = (101226.0-22764*s6)/5825.0 ;
  b16_11 = -169984.0/9087.0 ;
  b16_12 = -87.0/30290.0 ;
  b16_13 = 492.0/1165.0 ;
  b16_15 = 1260.0/233.0 ;
  e1 = -7.0/400.0 ;
  e8 = 63.0/200.0 ;
  e9 = -14.0/25.0 ;
  e10 = 21.0/20.0 ;
  e11 = -1024.0/975.0 ;
  e12 = -21.0/36400.0 ;
  e13 = -3.0/25.0 ;
  e14 = -9.0/280.0 ;
  e15 = 9.0/25.0 ;
  e16 = 233.0/4200.0 ;
  C36 = 1.0/36.0 ;
  om = 1.0 ;
  a = sqrt(2.0*om) ;
  ONES = arma::ones<vec>(x.n_rows, x.n_cols) ;
  u0 = a*(ONES/arma::cosh(sqrt(om)/g*x)) ;
  v = arma::fft(u0) ;
  tol = 1e-13 ;
  er1 = tol ;
  er2 = tol ;
  rho1 = 1.0 ;
  rho2 = 1.0 ;
  t = 0.0 ;
  Tf = 30.0 ;
  tloc = 0.0 ;
  tdraw = 0.1 ;
  dt2 = 1.0e-3 ;
  while (t<Tf)
  {
    er1 = 2*tol ;
    while (er1>tol)
    {
      rho1 = rho2 ;
      dt = dt2 ;
      k1 = RHS(v, 0, e, j, g, k, op, K2i) ;
      k2 = RHS(v+k1*a2*dt, a2*dt, e, j, g, k, op, K2i) ;
      k3 = RHS(v+(k1*b31+k2*b32)*dt, a3*dt, e, j, g, k, op, K2i) ;
      k4 = RHS(v+(k1*b41+k3*b43)*dt, a4*dt, e, j, g, k, op, K2i) ;
      k5 = RHS(v+(k1*b51+k3*b53+k4*b54)*dt, a5*dt, e, j, g, k, op, K2i) ;
      k6 = RHS(v+(k1*b61+k4*b64+k5*b65)*dt, a6*dt, e, j, g, k, op, K2i) ;
      k7 = RHS(v+(k1*b71+k4*b74+k5*b75+k6*b76)*dt, a7*dt, e, j, g, k, op, K2i) ;
      k8 = RHS(v+(k1*b81+k6*b86+k7*b87)*dt, a8*dt, e, j, g, k, op, K2i) ;
      k9 = RHS(v+(k1*b91+k6*b96+k7*b97+k8*b98)*dt, a9*dt, e, j, g, k, op, K2i) ;
      k10 = RHS(v+(k1*b10_1+k6*b10_6+k7*b10_7+k8*b10_8+k9*b10_9)*dt, a10*dt, e, j, g, k, op, K2i) ;
      k11 = RHS(v+(k1*b11_1+k7*b11_7+k8*b11_8+k9*b11_9+k10*b11_10)*dt, a11*dt, e, j, g, k, op, K2i) ;
      k12 = RHS(v+(k1*b12_1+k6*b12_6+k7*b12_7+k8*b12_8+k9*b12_9+k10*b12_10+k11*b12_11)*dt, a12*dt, e, j, g, k, op, K2i) ;
      k13 = RHS(v+(k1*b13_1+k6*b13_6+k7*b13_7+k8*b13_8+k9*b13_9+k10*b13_10+k11*b13_11+k12*b13_12)*dt, a13*dt, e, j, g, k, op, K2i) ;
      k14 = RHS(v+(k1*b14_1+k6*b14_6+k7*b14_7+k8*b14_8+k9*b14_9+k10*b14_10+k11*b14_11+k12*b14_12+k13*b14_13)*dt, dt, e, j, g, k, op, K2i) ;
      k15 = RHS(v+(k1*b15_1+k6*b15_6+k7*b15_7+k8*b15_8+k9*b15_9+k10*b15_10+k11*b15_11+k12*b15_12+k13*b15_13)*dt, a15*dt, e, j, g, k, op, K2i) ;
      k16 = RHS(v+(k1*b16_1+k6*b16_6+k7*b16_7+k8*b16_8+k9*b16_9+k10*b16_10+k11*b16_11+k12*b16_12+k13*b16_13+k15*b16_15)*dt, dt, e, j, g, k, op, K2i) ;
      er2 = dt*norm(k1*e1+k8*e8+k9*e9+k10*e10+k11*e11+k12*e12+k13*e13+k14*e14+k15*e15+k16*e16, "inf") ;
      rho2 = pow((tol/er2), C36)*pow((tol/er1), C36)*(pow(rho1, (-0.25))) ;
      dt2 = std::max((1+atan(rho2-1))*dt, 1e-6) ;
      er1 = er2 ;
    }
    w = v+(k1*c1+k8*c8+k9*c9+k10*c10+k11*c11+k12*c12+k13*c13+k14*c14)*dt ;
    v = turn(w, -dt, K2i, j) ;
    if ((tloc+dt>tdraw))
    {
      dt = tdraw-tloc ;
    }
    if ((t+dt>Tf))
    {
      dt = Tf-t ;
    }
    tloc = tloc+dt ;
    t = t+dt ;
    cout<<" t = "<<t<<endl;
    if (((tloc==tdraw) || (t==Tf)))
    {
      tloc = 0. ;
      u_hat = v ;
      u = arma::ifft<cx_vec>(u_hat) ;
    }
  }
  return 0 ;
}
