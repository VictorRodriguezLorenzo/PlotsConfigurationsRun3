#ifndef TOPDM_SPIN_OBSERVABLES_H
#define TOPDM_SPIN_OBSERVABLES_H

#include <cmath>
#include <limits>

#include "TLorentzVector.h"
#include "TVector3.h"
#include "ROOT/RVec.hxx"

using RVecF = ROOT::VecOps::RVec<float>;

namespace topdm {

enum SpinObservableIndex {
    kCosPhi = 0,
    kDevt = 1,
    kXkk = 2,
    kXrr = 3,
    kXnn = 4,
    kXkr = 5,
    kXkn = 6,
    kXrk = 7,
    kXrn = 8,
    kXnk = 9,
    kXnr = 10,
    kMttReco = 11,
    kPtttReco = 12,
    kYttReco = 13,
    kAbsDeltaYtt = 14
};

inline float spinNaN() {
    return std::numeric_limits<float>::quiet_NaN();
}

inline float clampUnit(double x) {
    if (x > 1.0) return 1.0f;
    if (x < -1.0) return -1.0f;
    return static_cast<float>(x);
}

/*
  Spin observables for a reconstructed dileptonic ttbar system.

  Required convention:
    - top     is the reconstructed t candidate associated with l+
    - antitop is the reconstructed tbar candidate associated with l-

  The helicity-frame axes follow the convention used in the paper:
      k = top direction in the ttbar rest frame
      n = (p x k) / |p x k|
      r = k x n
  with the pp beam direction p chosen so that p.k >= 0 (theta in [0, pi/2]).

  The charged leptons are first boosted to the reconstructed ttbar CM frame and
  then to their respective parent-top rest frames.  The second boosts are along
  +/-k, so the common k,r,n axes retain their spatial orientation.

  Output map:
     0 cosPhi
     1 Devt = -3*cosPhi
     2 xkk
     3 xrr
     4 xnn
     5 xkr
     6 xkn
     7 xrk
     8 xrn
     9 xnk
    10 xnr
    11 mtt_reco
    12 pttt_reco
    13 ytt_reco
    14 absDeltaY_tt
*/
inline RVecF spinObservables(
    const TLorentzVector& topLab,
    const TLorentzVector& antitopLab,
    const TLorentzVector& lplusLab,
    const TLorentzVector& lminusLab
) {
    RVecF out(15, spinNaN());

    const TLorentzVector ttLab = topLab + antitopLab;
    if (!std::isfinite(ttLab.E()) || ttLab.E() <= 0.0 ||
        !std::isfinite(ttLab.M2()) || ttLab.M2() <= 0.0) {
        return out;
    }

    out[kMttReco] = static_cast<float>(ttLab.M());
    out[kPtttReco] = static_cast<float>(ttLab.Pt());
    out[kYttReco] = static_cast<float>(ttLab.Rapidity());
    out[kAbsDeltaYtt] = static_cast<float>(
        std::abs(topLab.Rapidity() - antitopLab.Rapidity())
    );

    // 1) Boost everything to the reconstructed ttbar centre-of-mass frame.
    const TVector3 betaTT = -ttLab.BoostVector();
    if (!std::isfinite(betaTT.X()) || !std::isfinite(betaTT.Y()) ||
        !std::isfinite(betaTT.Z()) || betaTT.Mag2() >= 1.0) {
        return out;
    }

    TLorentzVector topCM = topLab;
    TLorentzVector antitopCM = antitopLab;
    TLorentzVector lpCM = lplusLab;
    TLorentzVector lmCM = lminusLab;

    topCM.Boost(betaTT);
    antitopCM.Boost(betaTT);
    lpCM.Boost(betaTT);
    lmCM.Boost(betaTT);

    if (topCM.Vect().Mag2() <= 0.0 || antitopCM.Vect().Mag2() <= 0.0)
        return out;

    const TVector3 k = topCM.Vect().Unit();

    // Boost the two proton beam directions into the ttbar CM.  The arbitrary
    // beam energy cancels when the spatial vector is normalized.
    TLorentzVector beamPlus(0.0, 0.0, +1.0, 1.0);
    TLorentzVector beamMinus(0.0, 0.0, -1.0, 1.0);
    beamPlus.Boost(betaTT);
    beamMinus.Boost(betaTT);

    if (beamPlus.Vect().Mag2() <= 0.0 || beamMinus.Vect().Mag2() <= 0.0)
        return out;

    const TVector3 pPlus = beamPlus.Vect().Unit();
    const TVector3 pMinus = beamMinus.Vect().Unit();
    const TVector3 p = (pPlus.Dot(k) >= pMinus.Dot(k)) ? pPlus : pMinus;

    TVector3 n = p.Cross(k);
    if (n.Mag2() < 1e-12) {
        // Stable orthogonal fallback for the measure-zero collinear case.
        TVector3 trial(1.0, 0.0, 0.0);
        if (std::abs(trial.Dot(k)) > 0.9)
            trial.SetXYZ(0.0, 1.0, 0.0);
        n = trial.Cross(k);
    }
    if (n.Mag2() <= 0.0) return out;
    n = n.Unit();

    TVector3 r = k.Cross(n);
    if (r.Mag2() <= 0.0) return out;
    r = r.Unit();

    // 2) From ttbar CM to each parent-top rest frame.
    const TVector3 betaTop = -topCM.BoostVector();
    const TVector3 betaAntiTop = -antitopCM.BoostVector();
    if (betaTop.Mag2() >= 1.0 || betaAntiTop.Mag2() >= 1.0)
        return out;

    TLorentzVector lpTop = lpCM;
    TLorentzVector lmAntiTop = lmCM;
    lpTop.Boost(betaTop);
    lmAntiTop.Boost(betaAntiTop);

    if (lpTop.Vect().Mag2() <= 0.0 || lmAntiTop.Vect().Mag2() <= 0.0)
        return out;

    const TVector3 up = lpTop.Vect().Unit();
    const TVector3 um = lmAntiTop.Vect().Unit();

    const float cpK = clampUnit(up.Dot(k));
    const float cpR = clampUnit(up.Dot(r));
    const float cpN = clampUnit(up.Dot(n));
    const float cmK = clampUnit(um.Dot(k));
    const float cmR = clampUnit(um.Dot(r));
    const float cmN = clampUnit(um.Dot(n));

    out[kXkk] = cpK * cmK;
    out[kXrr] = cpR * cmR;
    out[kXnn] = cpN * cmN;
    out[kXkr] = cpK * cmR;
    out[kXkn] = cpK * cmN;
    out[kXrk] = cpR * cmK;
    out[kXrn] = cpR * cmN;
    out[kXnk] = cpN * cmK;
    out[kXnr] = cpN * cmR;

    // By completeness of the orthonormal k,r,n basis this is the opening-angle
    // observable used in the paper.
    const float cosPhi = clampUnit(
        static_cast<double>(out[kXkk] + out[kXrr] + out[kXnn])
    );
    out[kCosPhi] = cosPhi;
    out[kDevt] = -3.0f * cosPhi;

    return out;
}

} // namespace topdm

#endif
