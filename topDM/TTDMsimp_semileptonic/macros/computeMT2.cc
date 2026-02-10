#ifndef COMPUTE_MT2_CC
#define COMPUTE_MT2_CC

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "TLorentzVector.h"

/*
 *  MT2 implementation based on the algorithm in:
 *  http://arxiv.org/abs/1411.4312
 *
 *  This code is adapted from the public implementation by Christopher Lester
 *  (University of Cambridge). See the original header for full copyright
 *  information and usage notes.
 */

namespace Lester {

struct EllipseParams {
  EllipseParams() = default;

  EllipseParams(const double c_xx2, const double c_yy2, const double c_xy2,
                const double c_x2, const double c_y2, const double c2)
      : c_xx(c_xx2), c_yy(c_yy2), c_xy(c_xy2), c_x(c_x2), c_y(c_y2), c(c2) {
    if (c_xx < 0 || c_yy < 0) {
      throw "precondition violation";
    }
    setDet();
  }

  EllipseParams(const double x0, const double y0)
      : c_xx(1),
        c_yy(1),
        c_xy(0),
        c_x(-x0),
        c_y(-y0),
        c(x0 * x0 + y0 * y0),
        det(0) {}

  void setDet() {
    det = (2.0 * c_x * c_xy * c_y + c * c_xx * c_yy - c_yy * c_x * c_x -
           c * c_xy * c_xy - c_xx * c_y * c_y);
  }

  double lesterFactor(const EllipseParams &e2) const {
    const EllipseParams &e1 = *this;
    const double ans = e1.c_xx * e1.c_yy * e2.c + 2.0 * e1.c_xy * e1.c_y * e2.c_x -
                       2.0 * e1.c_x * e1.c_yy * e2.c_x + e1.c * e1.c_yy * e2.c_xx -
                       2.0 * e1.c * e1.c_xy * e2.c_xy + 2.0 * e1.c_x * e1.c_y * e2.c_xy +
                       2.0 * e1.c_x * e1.c_xy * e2.c_y - 2.0 * e1.c_xx * e1.c_y * e2.c_y +
                       e1.c * e1.c_xx * e2.c_yy - e2.c_yy * (e1.c_x * e1.c_x) -
                       e2.c * (e1.c_xy * e1.c_xy) - e2.c_xx * (e1.c_y * e1.c_y);
    return ans;
  }

  bool operator==(const EllipseParams &other) const {
    return c_xx == other.c_xx && c_yy == other.c_yy && c_xy == other.c_xy &&
           c_x == other.c_x && c_y == other.c_y && c == other.c;
  }

  double c_xx = 0;
  double c_yy = 0;
  double c_xy = 0;
  double c_x = 0;
  double c_y = 0;
  double c = 0;
  double det = 0;  // determinant of 3x3 conic matrix
};

bool __private_ellipsesAreDisjoint(const double coeffLamPow3,
                                   const double coeffLamPow2,
                                   const double coeffLamPow1,
                                   const double coeffLamPow0) {
  if (fabs(coeffLamPow3) == 0) {
    throw 1;  // singular ellipse
  }

  const double a = coeffLamPow2 / coeffLamPow3;
  const double b = coeffLamPow1 / coeffLamPow3;
  const double c = coeffLamPow0 / coeffLamPow3;

  const double thing1 = -3.0 * b + a * a;
  if (thing1 <= 0) return false;
  const double thing2 = -27.0 * c * c + 18.0 * c * a * b + a * a * b * b -
                        4.0 * a * a * a * c - 4.0 * b * b * b;
  if (thing2 <= 0) return false;

  return ((a >= 0 && 3.0 * a * c + b * a * a - 4.0 * b * b < 0) ||
          (a < 0));
}

bool ellipsesAreDisjoint(const EllipseParams &e1, const EllipseParams &e2) {
  if (e1 == e2) {
    return false;
  }

  const double coeffLamPow3 = e1.det;
  const double coeffLamPow2 = e1.lesterFactor(e2);
  const double coeffLamPow1 = e2.lesterFactor(e1);
  const double coeffLamPow0 = e2.det;

  if (fabs(coeffLamPow3) >= fabs(coeffLamPow0)) {
    return __private_ellipsesAreDisjoint(coeffLamPow3, coeffLamPow2,
                                         coeffLamPow1, coeffLamPow0);
  }
  return __private_ellipsesAreDisjoint(coeffLamPow0, coeffLamPow1,
                                       coeffLamPow2, coeffLamPow3);
}

}  // namespace Lester

class asymm_mt2_lester_bisect {
 public:
  static const int MT2_ERROR = -1;

  static double get_mT2(const double mVis1, const double pxVis1,
                        const double pyVis1, const double mVis2,
                        const double pxVis2, const double pyVis2,
                        const double pxMiss, const double pyMiss,
                        const double mInvis1, const double mInvis2,
                        const double desiredPrecisionOnMT2 = 0,
                        const bool useDeciSectionsInitially = true) {
    const double mT2_Sq = get_mT2_Sq(
        mVis1, pxVis1, pyVis1, mVis2, pxVis2, pyVis2, pxMiss, pyMiss,
        mInvis1, mInvis2, desiredPrecisionOnMT2, useDeciSectionsInitially);
    if (mT2_Sq == MT2_ERROR) {
      return MT2_ERROR;
    }
    return sqrt(mT2_Sq);
  }

  static double get_mT2_Sq(const double mVis1, const double pxVis1,
                           const double pyVis1, const double mVis2,
                           const double pxVis2, const double pyVis2,
                           const double pxMiss, const double pyMiss,
                           const double mInvis1, const double mInvis2,
                           const double desiredPrecisionOnMT2 = 0,
                           const bool useDeciSectionsInitially = true) {
    const double m1Min = mVis1 + mInvis1;
    const double m2Min = mVis2 + mInvis2;

    if (m1Min > m2Min) {
      return asymm_mt2_lester_bisect::get_mT2_Sq(
          mVis2, pxVis2, pyVis2, mVis1, pxVis1, pyVis1, pxMiss, pyMiss,
          mInvis2, mInvis1, desiredPrecisionOnMT2, useDeciSectionsInitially);
    }

    assert(m1Min <= m2Min);
    const double mMin = m2Min;

    const double msSq = mVis1 * mVis1;
    const double sx = pxVis1;
    const double sy = pyVis1;
    const double mpSq = mInvis1 * mInvis1;

    const double mtSq = mVis2 * mVis2;
    const double tx = pxVis2;
    const double ty = pyVis2;
    const double mqSq = mInvis2 * mInvis2;

    const double sSq = sx * sx + sy * sy;
    const double tSq = tx * tx + ty * ty;
    const double pMissSq = pxMiss * pxMiss + pyMiss * pyMiss;
    const double massSqSum = msSq + mtSq + mpSq + mqSq;
    const double scaleSq = (massSqSum + sSq + tSq + pMissSq) / 8.0;

    if (scaleSq == 0) {
      return 0;
    }
    const double scale = sqrt(scaleSq);

    double mLower = mMin;
    double mUpper = mMin + scale;
    unsigned int attempts = 0;
    const unsigned int maxAttempts = 10000;
    while (true) {
      ++attempts;
      const double mUpperSq = mUpper * mUpper;
      const Lester::EllipseParams &side1 = helper(mUpperSq, msSq, -sx, -sy,
                                                  mpSq, 0, 0);
      const Lester::EllipseParams &side2 = helper(mUpperSq, mtSq, +tx, +ty,
                                                  mqSq, pxMiss, pyMiss);

      bool disjoint;
      try {
        disjoint = Lester::ellipsesAreDisjoint(side1, side2);
      } catch (...) {
        return MT2_ERROR;
      }

      if (!disjoint) {
        break;
      }

      if (attempts >= maxAttempts) {
        std::cerr << "MT2 algorithm failed to find upper bound to MT2"
                  << std::endl;
        return MT2_ERROR;
      }

      mUpper *= 2;
    }

    bool goLow = useDeciSectionsInitially;
    while (desiredPrecisionOnMT2 <= 0 || mUpper - mLower > desiredPrecisionOnMT2) {
      const double trialM = (goLow ? (mLower * 15 + mUpper) / 16
                                   : (mUpper + mLower) / 2.0);

      if (trialM <= mLower || trialM >= mUpper) {
        return trialM * trialM;
      }

      const double trialMSq = trialM * trialM;
      const Lester::EllipseParams &side1 = helper(trialMSq, msSq, -sx, -sy,
                                                  mpSq, 0, 0);
      const Lester::EllipseParams &side2 = helper(trialMSq, mtSq, +tx, +ty,
                                                  mqSq, pxMiss, pyMiss);

      try {
        const bool disjoint = Lester::ellipsesAreDisjoint(side1, side2);
        if (disjoint) {
          mLower = trialM;
          goLow = false;
        } else {
          mUpper = trialM;
        }
      } catch (...) {
        return mLower * mLower;
      }
    }

    const double mAns = (mLower + mUpper) / 2.0;
    return mAns * mAns;
  }

 private:
  static const Lester::EllipseParams helper(const double mSq,
                                            const double mtSq, const double tx,
                                            const double ty, const double mqSq,
                                            const double pxmiss,
                                            const double pymiss) {
    const double txSq = tx * tx;
    const double tySq = ty * ty;
    const double pxmissSq = pxmiss * pxmiss;
    const double pymissSq = pymiss * pymiss;

    const double c_xx = +4.0 * mtSq + 4.0 * tySq;
    const double c_yy = +4.0 * mtSq + 4.0 * txSq;
    const double c_xy = -4.0 * tx * ty;
    const double c_x = -4.0 * mtSq * pxmiss - 2.0 * mqSq * tx +
                       2.0 * mSq * tx - 2.0 * mtSq * tx + 4.0 * pymiss * tx * ty -
                       4.0 * pxmiss * tySq;

    const double c_y = -4.0 * mtSq * pymiss - 4.0 * pymiss * txSq -
                       2.0 * mqSq * ty + 2.0 * mSq * ty - 2.0 * mtSq * ty +
                       4.0 * pxmiss * tx * ty;

    const double c = -mqSq * mqSq + 2 * mqSq * mSq - mSq * mSq + 2 * mqSq * mtSq +
                     2 * mSq * mtSq - mtSq * mtSq + 4.0 * mtSq * pxmissSq +
                     4.0 * mtSq * pymissSq + 4.0 * mqSq * pxmiss * tx -
                     4.0 * mSq * pxmiss * tx + 4.0 * mtSq * pxmiss * tx +
                     4.0 * mqSq * txSq + 4.0 * pymissSq * txSq +
                     4.0 * mqSq * pymiss * ty - 4.0 * mSq * pymiss * ty +
                     4.0 * mtSq * pymiss * ty - 8.0 * pxmiss * pymiss * tx * ty +
                     4.0 * mqSq * tySq + 4.0 * pxmissSq * tySq;

    return Lester::EllipseParams(c_xx, c_yy, c_xy, c_x, c_y, c);
  }
};

// ---------------------------------------------------------------------------
// Compute MT2 from two leptons and MET using the Lester bisect algorithm
// ---------------------------------------------------------------------------
float computeMT2(float l1_pt, float l1_eta, float l1_phi, float l2_pt,
                 float l2_eta, float l2_phi, float met_pt, float met_phi) {
  TLorentzVector lepton1, lepton2;
  lepton1.SetPtEtaPhiM(l1_pt, l1_eta, l1_phi, 0.0);
  lepton2.SetPtEtaPhiM(l2_pt, l2_eta, l2_phi, 0.0);

  const double mVisA = fabs(lepton1.M());
  const double mVisB = fabs(lepton2.M());

  const double pxA = lepton1.Px();
  const double pyA = lepton1.Py();
  const double pxB = lepton2.Px();
  const double pyB = lepton2.Py();

  const double pxMiss = met_pt * cos(met_phi);
  const double pyMiss = met_pt * sin(met_phi);

  const double chiA = 0.0;
  const double chiB = 0.0;
  const double desiredPrecisionOnMt2 = 0.0;

  const double mT2 = asymm_mt2_lester_bisect::get_mT2(
      mVisA, pxA, pyA, mVisB, pxB, pyB, pxMiss, pyMiss, chiA, chiB,
      desiredPrecisionOnMt2);

  return static_cast<float>(mT2);
}

#endif

