#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <arm_neon.h>

float rumach() {
  float comp = 0.0e0;
  float u = 1.0e0;

  while (comp != 1.0e0) {
    u = u*0.5e0;
    comp = 1.0e0 + u;
  }

  return u*2.0e0;
}

/*
 * SINTDY
 */

void sintdy(
    float t, int n, int k,
    float yh[4][2], float* dky,
    float h, float hu, float tn,
    int nyh, int iflag,
    int nq, int l
) {
  float uround = rumach();
  iflag = 0;

  if ((k < 0) | (k > nq)) {
    std::cout << "Error: SINTDY- K (=I1) illegal";
    iflag = -1;
    return;
  }

  float tp = tn - hu - 100.0e0*uround*copysign(abs(tn) + abs(hu), hu);
  if ((t-tp)*(t-tn) > 0.0e0) {
    std::cout << "Error: SINTDY- T (=R1) illegal";
    iflag = -2;
    return;
  }

  float s = (t - tn)/h;
  int ic = 1;
  float c;
  
  if (k != 0) {
    int jj1 = l - k;

    for (int jj = jj1; jj < nq; jj++) {
      ic = ic * jj;
    }
  }

  c = ic;

  for (int i = 0; i < n; i++) {
    dky[i] = c*yh[i][l];
  }

  if (k == nq) {
    float r = pow(h, -k);
    for (int i = 0; i < n; i+=4) {
      float32x4_t dky_v = vld1q_f32(dky+i);
      vst1q_f32(dky+i, vmulq_n_f32(dky_v, r));
    }

    return;
  }

  int jb2 = nq - k;
  for (int jb = 1; jb < jb2; jb++) {
    int j = nq - jb;
    int jp1 = j + 1;
    ic = 1;

    if (k != 0) {
      int jj1 = jp1 - k;

      for (int jj = jj1; jj < j; jj++) {
        ic = ic*jj;
      }
    }

    c = ic;

    for (int i = 0; i < n; i+=4) {
      float32x4_t dky_v = vld1q_f32(dky+i);
      vst1q_f32(dky+i, vmulq_n_f32(dky_v, s));

      dky[i] += c*yh[i][jp1];
      dky[i+1] += c*yh[i+1][jp1];
      dky[i+2] += c*yh[i+2][jp1];
      dky[i+3] += c*yh[i+3][jp1];
    }
  }

  if (k == 0) { return; }
  
  float r = pow(h, -k);
  for (int i = 0; i < n; i+=4) {
    float32x4_t dky_v = vld1q_f32(dky+i);
    vst1q_f32(dky+i, vmulq_n_f32(dky_v, r));
  }

  return;
}

int main() {
  // initialize arguments
  int k = 2;
  int n = 4;
  int l = 2;
  int nq = 4;
  int nyh = n;
  int iflag = 0;
  float t = 0.4;
  float tn = t;
  float h = 1.0e0;
  float hu = 0.0e0;
  float dky[n] = {1.0, 0.0, 0.0, 0.0};
  float yh[4][2];
  yh[0][0] = 1.0;
  yh[0][1] = 0.0;
  yh[1][0] = 0.0;
  yh[1][1] = 1.0;
  yh[2][0] = 1.0;
  yh[2][1] = 1.0;
  yh[3][0] = 0.0;
  yh[3][1] = 0.0;

  sintdy(t, n, k, yh, dky, h, hu, tn, nyh, iflag, nq, l);

  // print dky
  std::cout << "\nDKY[0] = " << dky[0]
    << "\nDKY[1] = " << dky[1]
    << "\nDKY[2] = " << dky[2]
    << "\nDKY[3] = " << dky[3];

  return 0;
}
