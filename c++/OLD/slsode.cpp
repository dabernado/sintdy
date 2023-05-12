#include <iostream>
#include <functional>
#include <stdlib>
#include <math>

void main() {
  // initialize arguments
  int neq = 3;
  int t = 0;
  int itol = 2;
  int itask = 1;
  int istate = 1;
  int iopt = 0;
  int lrw = 58;
  int liw = 23;
  int mf = 21;
  int iwork[23] = 0;
  float t_out = 0.4;
  float rtol = 1.e-4;
  float y[3] = {1.0, 0.0, 0.0};
  float atol[3] = {1.e-6, 1.e-10, 1.e-6};

  for (int i = 0; i < 12; i++) {
    // call function
    slsode(
        fex, jex,
        neq, itol, itask, istate, lrw, liw, mf,
        iwork,
        y, atol,
        t, t_out, rtol,
        iopt);

    // print results
    cout << "\nAt t = " << t
      << "  y = " << y[0] << " " << y[1] << " " << y[2];

    // check istate for error
    if (istate < 0) {
      cout << "\nError halt.. ISTATE = " << istate;
      return;
    }

    t_out = t_out * 10;
  }

  // print iwork
  cout << "\nNo. steps = " << iwork[10]
    << " No. f-s = " << iwork[11]
    << " No. J-s = " << iwork[12];

  return;
}

void fex(int neq, float t, float y[3], float ydot[3]) {
  ydot[0] = (-0.04 * y[0]) + (1.e4 * y[1] * y[2]);
  ydot[2] = 3.e7 * y[1] * y[1];
  ydot[1] = -(ydot[0]) - ydot[2];
  return;
}

void jex(
    int neq,
    float t,
    float y[3],
    int ml, int mu, int nrpd,
    float pd[nrpd][3]
) {
  pd[0][0] = -0.04;
  pd[0][1] = 1.e4 * y[2];
  pd[0][2] = 1.e4 * y[1];
  pd[1][0] = 0.04;
  pd[1][2] = -(pd[0][2]);
  pd[2][1] = 6.e7 * y[1];
  pd[1][1] = -(pd[0][1]) - pd[2][1];
  return;
}

/*
 * SLSODE
 */

// add common variables as statics
void slsode(
    std::function<
      void(int, float, float*, float*)
      > f,
    std::function<
      void(int, int, int, int, float, float*, float*)
      > jac,

    int neq, int itol, int itask, int istate, int lrw, int liw, int mf,
    int iwork[23],
    float y[3], float atol[3],
    float t, float t_out, float rtol,
    bool iopt
) {
  static int maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu;
  static int icf, ierpj, iersl, jcur, jstart, kflag, l;
  static int lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter;
  static float ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround;
  static int init, mxstep, mxhnil, nhnil, nslast, nyh;
  static int iownd[6], iowns[6];
  static float rowns[209];

  const int mord[2] = {12, 5};
  const int mxstp0 = 500;
  const int mxhnl0 = 10;

  int h0, i1, i2;
  float tcrit;

  n = neq;
  meth = mf/10;
  miter = mf - (10*meth);
  int ml = iwork[0];
  int mu = iwork[1];

  // check if meth is within bounds of mord
  if ((meth < 1) | (meth > 2)) {
    cout << "Error: METH out of bounds";
    return;
  }

  if (iopt != 1) {
    maxord = mord[meth];
    mxstep = mxstp0;
    mxhnil = mxhnl0;

    if (istate == 1) {
      h0 = 0.0e0;
    }

    hmxi = 0.0e0;
    hmin = 0.0e0;
  } else {
    maxord = iwork[4];
    if (maxord == 0) { maxord = 100; }
    maxord = min(maxord, mord[meth]);
    
    mxstep = iwork[5];
    if (mxstep == 0) { mxstep = mxstp0; }
    
    mxhnil = iwork[6];
    if (mxhnil == 0) { mxhnil = mxhnl0; }

    if (istate == 1) {
      h0 = rwork[4];
    }

    hmax = rwork[5]
    hmxi = 0.0e0;
    if (hmax > 0.0e0) { hmxi = 1.0e0/hmax; }

    hmin = rwork[6];
  }

  lyh = 21;
  if (istate == 1) { nyh = n; }

  lwm = lyh + (maxord + 1)*nyh;

  if (miter == 0) {
    lenwm = 0;
  } else if ((miter == 1) | (miter == 2)) {
    lenwm = n*n + 2;
  } else if (miter == 3) {
    lenwm = n + 2;
  } else if (miter >= 4) {
    lenwm = (2*ml + mu + 1)*n + 2;
  }

  lewt = lwm + lenwm;
  lsavf = lewt + n;
  lacor = lsavf + n;
  lenrw = lacor + n - 1;
  iwork[16] = lenrw;
  liwm = 1;
  leniw = 20 + n;
  if ((miter == 0) | (miter == 3)) { leniw = 20; }
  iwork[17] = leniw;

  if (istate != 1) {
    jstart = -1;

    if (nq > maxord) {
      for (int i = 0; i < n; i++) {
        rwork[i+lsavf-1] = rwork[i+lwm-1];
      }
    }

    if (miter > 0) {
      rwork[lwm] = sqrt(uround);
    }

    if (n != nyh) {
      i1 = lyh + l*nyh;
      i2 = lyh + (maxord + 1)*nyh - 1;

      if (i1 <= i2) {
        for (int i = i1; i < i2; i++) {
          rwork[i] = 0.0e0;
        }
      }
    }

    nslast = nst;
    switch(itask) {
      case 0:
        if ((tn - t_out)*h >= 0.0e0) {
          sintdy(t_out, 0, rwork[lyh], nyh, y, iflag);
          t = t_out;
          goto L420;
        }

        goto default;
      case 2:
        tp = tn - hu*(1.0e0 + 100.0e0*uround);
        if ((tn - t_out)*h >= 0.0e0) {
          goto L400;
        }

        goto default;
      case 3:
        tcrit = rwork[0];
        if ((tn - t_out)*h < 0.0e0) {
          goto L245;
        }

        sintdy(t_out, 0, rwork[lyh], nyh, y, iflag);
        t = t_out;
        goto L420;
      case 4:
        tcrit = rwork[0];
        hmx = abs(tn) + abs(h);
L245:
        ihit = abs(tn - tcrit) <= 100.0e0*uround*hmx;
        if (ihit) { goto L400; }

        tnext = tn + h*(1.0e0 + 4.0e0*uround);
        if ((tnext - tcrit)*h <= 0.0e0) {
          goto default;
        }

        h = (tcrit - tn)*(1.0e0 - 4.0e0*uround);
        if (istate == 2) {
          jstart = -2;
        }
      default:
        break;
    }

L400:
L420:

  } else {
    uround = 1.4e-45; // set uround to smallest positive float
    tn = t;

    if ((itask == 4) | (itask == 5)) {
      tcrit = rwork[0];

      if ((h0 != 0.0e0) & (((t + h0 - tcrit) * h0) > 0.0e0)) {
        h0 = tcrit - t;
      }
    }

    jstart = 0;
    if (miter > 0) { rwork[lwm] = sqrt(uround); }
    nhnil = 0;
    nst = 0;
    nje = 0;
    nslast = 0;
    hu = 0.0e0;
    nqu = 0;
    ccmax = 0.3e0;
    maxcor = 3;
    msbp = 20;
    mxncf = 10;

    lf0 = lyh + nyh;
    f(neq, t, y, rwork[lf0]);
    nfe = 1;

    for (int i = 0; i < n; i++) {
      rwork[i+lyh-1] = y[i];
    }

    nq = 1;
    h = 1.0e0;
    sewset(n, itol, rtol, atol, rwork[lyh], rwork[lewt]);
    for (int i = 0; i < n; i++) {
      rwork[i+lewt-1] = 1.0e0/rwork[i+lewt-1];
    }

    if (h0 == 0.0e0) {
      tdist = abs(tout - t);
      w0 = max(abs(t), abs(tout));
      tol = rtol[0];

      if (itol > 2) {
        for (int i = 0; i < n; i++) {
          tol = max(tol, rtol[i]);
        }
      }

      if (tol <= 0.0e0) {
        atoli = atol[0];

        for (int i = 0; i < n; i++) {
          if ((itol == 2) | (itol == 4)) { atoli = atol[i]; }
          ayi = abs(y[i]);
          if (ayi != 0.0e0) { tol = max(tol, atoli/ayi); }
        }
      }

      tol = max(tol, 100.0e0*uround);
      tol = min(tol, 0.001e0);
      sum = svnorm(n, rwork[lf0], rwork[lewt]);
      sum = 1.0e0/(tol*w0*w0) + tol*pow(sum, 2);
      h0 = 1.0e0/sqrt(sum);
      h0 = min(h0, tdist);
      h0 = copysign(h0, tout - t);
    }

    rh = abs(h0)*hmxi;
    if (rh > 1.0e0) { h0 = h0/rh; }
    h = h0
    for (int i = 0; i < n; i++) {
      rwork[i+lf0-1] = h0*rwork[i+lf0-1];
    }
  }
}

void sintdy(float t, int k, float* yh, int nyh, float dky, int iflag) {
  iflag = 0;

  if ((k < 0) | (k > nq)) {
    cout << "Error: SINTDY- K (=I1) illegal";
    iflag = -1;
    return;
  }

  float tp = tn - hu - (100.0e0 * uround * copysign(abs(tn) + abs(hu), hu));
  if ((t-tp)*(t-tn) > 0.0e0) {
    cout << "Error: SINTDY- T (=R1) illegal";
    iflag = -2;
    return;
  }

  float s = (t - tn)/h;
  float c;
  int ic = 1;
  
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
    for (int i = 0; i < n; i++) {
      dky[i] = r*dky[i];
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

    for (int i = 0; i < n; i++) {
      dky[i] = (c * yh[i][jp1]) + (s * dky[i]);
    }
  }

  if (k == 0) { return; }
  
  float r = pow(h, -k);
  for (int i = 0; i < n; i++) {
    dky[i] = r * dky[i];
  }

  return;
}
