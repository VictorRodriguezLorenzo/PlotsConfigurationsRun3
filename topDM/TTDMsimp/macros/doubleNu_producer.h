#ifndef MKSHAPESRDF_PROCESSOR_MODULES_DOUBLENEUTRINOSOLUTIONMACRO_H
#define MKSHAPESRDF_PROCESSOR_MODULES_DOUBLENEUTRINOSOLUTIONMACRO_H

#include <Math/Factory.h>
#include <Math/IFunction.h>
#include <Math/Minimizer.h>
#include <Math/VectorUtil.h>
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector2.h>
#include <TVectorD.h>
#include <TVector3.h>
#include <TMatrixDSymEigen.h>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <memory>
#include <functional>
#include <string>
#include <cstdlib>

namespace nuana {

// ---------- Utilities ----------

// Toggle to force-enable the MET minimizer fallback when the ellipse
// intersection search fails.
constexpr bool kEnableMinimizerFallback = true;

// UnitCircle: returns a 3x3 matrix representing the unit circle in the F' coordinate system
inline TMatrixD UnitCircle() {
    TMatrixD U(3,3);
    U.Zero();
    U(0,0) = 1.0;
    U(1,1) = 1.0;
    U(2,2) = -1.0;
    return U;
}

// Construct a robust 2x2 MET covariance matrix
inline TMatrixD makeMetCov(double cov_xx, double cov_xy, double cov_yy) {
    TMatrixD sigma2(2, 2);
    sigma2(0, 0) = cov_xx;
    sigma2(0, 1) = cov_xy;
    sigma2(1, 0) = cov_xy;
    sigma2(1, 1) = cov_yy;

    bool finite = std::isfinite(cov_xx) && std::isfinite(cov_xy) && std::isfinite(cov_yy);
    double det = sigma2(0, 0) * sigma2(1, 1) - sigma2(0, 1) * sigma2(1, 0);
    double minDiag = std::min(sigma2(0, 0), sigma2(1, 1));

    if (!finite || minDiag <= 0.0 || det <= 0.0) {
        sigma2.UnitMatrix();
    }

    return sigma2;
}

inline double det3x3(const TMatrixD &M) {
    if (M.GetNrows() != 3 || M.GetNcols() != 3) {
        return 0.0;
    }
    const double a = M(0,0);
    const double b = M(0,1);
    const double c = M(0,2);
    const double d = M(1,0);
    const double e = M(1,1);
    const double f = M(1,2);
    const double g = M(2,0);
    const double h = M(2,1);
    const double i = M(2,2);
    return a * (e * i - f * h)
         - b * (d * i - f * g)
         + c * (d * h - e * g);
}

inline bool invert3x3(const TMatrixD &M, TMatrixD &inv, double tol = 1e-12) {

    double det = det3x3(M);
    if (!std::isfinite(det) || std::abs(det) <= tol) {
        return false;
    }

    inv.ResizeTo(3, 3);

    inv(0,0) =  (M(1,1) * M(2,2) - M(1,2) * M(2,1)) / det;
    inv(0,1) = -(M(0,1) * M(2,2) - M(0,2) * M(2,1)) / det;
    inv(0,2) =  (M(0,1) * M(1,2) - M(0,2) * M(1,1)) / det;
    inv(1,0) = -(M(1,0) * M(2,2) - M(1,2) * M(2,0)) / det;
    inv(1,1) =  (M(0,0) * M(2,2) - M(0,2) * M(2,0)) / det;
    inv(1,2) = -(M(0,0) * M(1,2) - M(0,2) * M(1,0)) / det;
    inv(2,0) =  (M(1,0) * M(2,1) - M(1,1) * M(2,0)) / det;
    inv(2,1) = -(M(0,0) * M(2,1) - M(0,1) * M(2,0)) / det;
    inv(2,2) =  (M(0,0) * M(1,1) - M(0,1) * M(1,0)) / det;


    return true;
}

inline bool invert2x2(const TMatrixD &M, TMatrixD &inv, double tol = 1e-12) {
    if (M.GetNrows() != 2 || M.GetNcols() != 2) {
        return false;
    }

    double det = M(0,0) * M(1,1) - M(0,1) * M(1,0);
    if (!std::isfinite(det) || std::abs(det) <= tol) {
        return false;
    }

    inv.ResizeTo(2, 2);
    const double invDet = 1.0 / det;
    inv(0,0) =  M(1,1) * invDet;
    inv(0,1) = -M(0,1) * invDet;
    inv(1,0) = -M(1,0) * invDet;
    inv(1,1) =  M(0,0) * invDet;
    return true;
}

// Cofactor for 3x3 matrix A at (i,j)
inline double cofactor(const TMatrixD &A, int i, int j) {
    // Compute determinant of minor (2x2) by skipping row i and column j
    int rows[2], cols[2], r = 0, c = 0;
    for (int idx = 0; idx < 3; ++idx) {
        if (idx != i) rows[r++] = idx;
        if (idx != j) cols[c++] = idx;
    }
    double det2 = A(rows[0], cols[0]) * A(rows[1], cols[1]) - A(rows[0], cols[1]) * A(rows[1], cols[0]);
    return ((i + j) % 2 == 0 ? 1.0 : -1.0) * det2;
}

inline TMatrixD Rotation(int axis, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);

    TMatrixD R(3, 3);
    R.UnitMatrix();
    R *= c;

    for (int i = -1; i <= 1; ++i) {
        const int row = (axis - i + 3) % 3;
        const int col = (axis + i + 3) % 3;
        R(row, col) = i * s + (1 - i * i);
    }

    return R;
}

inline TMatrixD Derivative() {
    // Matrix to differentiate [cos(theta), sin(theta), 1]
    double angle = M_PI / 2.0;
    double c = std::cos(angle);
    double s = std::sin(angle);
    TMatrixD rot(3,3);
    rot.UnitMatrix();
    rot(0,0) = c; rot(0,1) = -s;
    rot(1,0) = s; rot(1,1) = c;
    rot(2,2) = 1.0;

    TMatrixD diag(3,3);
    diag.Zero();
    diag(0,0) = 1.0;
    diag(1,1) = 1.0;
    diag(2,2) = 0.0;

    TMatrixD result = rot * diag;
    return result;
}

// Valid real sqrt solutions of y = x^2
inline std::vector<double> multisqrt(double y) {
    if (y < 0.0) {
        return {};
    } else if (y == 0.0) {
        return {0.0};
    } else {
        double r = std::sqrt(y);
        return {-r, r};
    }
}

// Check whether a homogeneous coordinate lies on a conic section defined by "conic".
inline bool satisfiesConic(const TVectorD &v, const TMatrixD &conic, double tol = 1e-6) {
    if (v.GetNrows() != 3 || conic.GetNrows() != 3 || conic.GetNcols() != 3) {
        return false;
    }

    TVectorD Cv = conic * v;
    double val = 0.0;
    double norm = 0.0;

    for (int i = 0; i < 3; ++i) {
        val += v(i) * Cv(i);
        norm += v(i) * v(i);
    }

    if (!std::isfinite(val) || !std::isfinite(norm)) {
        return false;
    }

    double scale = std::max(norm, 1.0);
    double normalized = val / scale;
    return std::abs(normalized) <= tol;
}

// factor_degenerate: linear factors (lines) for degenerate quadratic (3x3 symmetric G)
inline std::vector<std::array<double,3>> factor_degenerate(
    const TMatrixD &G, double zero = 0.0
) {
    std::vector<std::array<double,3>> lines;

    if (std::abs(G(0,0)) <= zero && std::abs(G(1,1)) <= zero) {
        lines.push_back({G(0,1), 0.0, G(1,2)});
        lines.push_back({0.0, G(0,1), G(0,2) - G(1,2)});
        return lines;
    }

    bool swapXY = std::abs(G(0,0)) > std::abs(G(1,1));
    TMatrixD Q(G);
    if (swapXY) {
        TMatrixD tmp(3,3);
        int order[3] = {1, 0, 2};
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                tmp(r,c) = Q(order[r], order[c]);
        Q = tmp;
    }

    double denom = Q(1,1);
    if (std::abs(denom) <= zero) {
        return lines;
    }

    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            Q(r,c) /= denom;

    double q22 = cofactor(Q, 2, 2);

    auto swap_back = [](const std::array<double,3> &L) {
        return std::array<double,3>{L[1], L[0], L[2]};
    };

    if (-q22 <= zero) {
        double cof00 = cofactor(Q, 0, 0);
        auto roots = multisqrt(-cof00);
        for (double sVal : roots) {
            std::array<double,3> L{Q(0,1), Q(1,1), Q(1,2) + sVal};
            if (swapXY) {
                L = swap_back(L);
            }
            lines.push_back(L);
        }
    } else {
        double x0 = cofactor(Q, 0, 2) / q22;
        double y0 = cofactor(Q, 1, 2) / q22;
        auto roots = multisqrt(-q22);
        for (double sVal : roots) {
            double m = Q(0,1) + sVal;
            std::array<double,3> L{m, Q(1,1), -Q(1,1) * y0 - m * x0};
            if (swapXY) {
                L = swap_back(L);
            }
            lines.push_back(L);
        }
    }

    return lines;
}

inline std::vector<TVectorD> intersections_ellipse_line(
    const TMatrixD &ellipse,
    const std::array<double,3> &line,
    double zero = 1e-12
) {
    TMatrixD cross(3,3);
    for (int i = 0; i < 3; ++i) {
        TVectorD row(3);
        for (int j = 0; j < 3; ++j) {
            row[j] = ellipse(i,j);
        }
        cross(i,0) = line[1]*row[2] - line[2]*row[1];
        cross(i,1) = line[2]*row[0] - line[0]*row[2];
        cross(i,2) = line[0]*row[1] - line[1]*row[0];
    }

    TMatrixD crossT(TMatrixD::kTransposed, cross);
    TMatrixDEigen eig(crossT);
    const TMatrixD &eigVecs = eig.GetEigenVectors();

    std::vector<std::pair<TVectorD,double>> candidates;

    for (int i = 0; i < 3; ++i) {
        TVectorD v(3);
        for (int j = 0; j < 3; ++j) {
            v[j] = eigVecs(j, i);
        }

        if (std::abs(v[2]) <= zero) {
            continue;
        }

        double inv = 1.0 / v[2];
        TVectorD s_vec(3);
        s_vec[0] = v[0] * inv;
        s_vec[1] = v[1] * inv;
        s_vec[2] = 1.0;

        double lv = line[0]*v[0] + line[1]*v[1] + line[2]*v[2];
        TVectorD Ev = ellipse * v;
        double vev = v[0]*Ev[0] + v[1]*Ev[1] + v[2]*Ev[2];

        double k = lv*lv + vev*vev;
        candidates.emplace_back(s_vec, k);
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const auto &a, const auto &b) { return a.second < b.second; });

    if (candidates.size() > 2) {
        candidates.resize(2);
    }

    std::vector<TVectorD> result;
    for (const auto &entry : candidates) {
        if (entry.second < zero) {
            result.push_back(entry.first);
        }
    }
    return result;
}

inline std::pair<std::vector<TVectorD>, std::vector<std::array<double,3>>>
intersections_ellipses(
    const TMatrixD &A, const TMatrixD &B,
    bool returnLines = false, double zero = 1e-10
) {
    double detA = det3x3(A);
    double detB = det3x3(B);

    const TMatrixD *AA = &A;
    const TMatrixD *BB = &B;
    if (std::abs(detB) > std::abs(detA)) {
        AA = &B;
        BB = &A;
    }

    TMatrixD invA;
    if (!invert3x3(*AA, invA)) {
        return {std::vector<TVectorD>{}, std::vector<std::array<double,3>>{}};
    }

    TMatrixD M = invA * (*BB);
    TMatrixDEigen eig(M);
    TVectorD evalsRe = eig.GetEigenValuesRe();
    TVectorD evalsIm = eig.GetEigenValuesIm();

    double eigenvalue = 0.0;
    bool found = false;
    for (int i = 0; i < evalsRe.GetNrows(); ++i) {
        if (std::abs(evalsIm(i)) <= zero) {
            eigenvalue = evalsRe(i);
            found = true;
            break;
        }
    }
    if (!found) {
        return {std::vector<TVectorD>{}, std::vector<std::array<double,3>>{}};
    }

    TMatrixD G = (*BB) - eigenvalue * (*AA);
    auto lines = factor_degenerate(G, zero);

    std::vector<TVectorD> points;
    for (const auto &line : lines) {
        auto pts = intersections_ellipse_line(*AA, line, zero);
        points.insert(points.end(), pts.begin(), pts.end());
    }

    if (returnLines) {
        return {points, lines};
    }
    return {points, std::vector<std::array<double,3>>{}};
}

struct nuSolutionSet {
    TLorentzVector b, mu;
    double c, s, x0, x0p, Sx, Sy, w, w_, x1, y1, Z, Om2, eps2, mW2;
    mutable bool usedMatrixFallback_;

    nuSolutionSet()
        : b(), mu(), c(1.0), s(0.0), x0(0.0), x0p(0.0), Sx(0.0), Sy(0.0), w(0.0), w_(0.0),
          x1(0.0), y1(0.0), Z(0.0), Om2(1.0), eps2(0.0), mW2(80.385 * 80.385),
          usedMatrixFallback_(false)
    {
        b.SetPxPyPzE(0.0, 0.0, 0.0, 1.0);
        mu.SetPxPyPzE(0.0, 0.0, 0.0, 1.0);
    }

    nuSolutionSet(const TLorentzVector& b_, const TLorentzVector& mu_,
          double mW = 80.385, double mT = 172.5, double mN = 0.0)
        : b(b_), mu(mu_), usedMatrixFallback_(false)
    {
        mW2 = mW * mW;
        double mT2 = mT * mT;
        double mN2 = mN * mN;

        c = ROOT::Math::VectorUtil::CosTheta(b, mu);
        s = std::sqrt(1.0 - c * c);

        x0p = - (mT2 - mW2 - b.M2()) / (2.0 * b.E());
        x0  = - (mW2 - mu.M2() - mN2) / (2.0 * mu.E());

        double Bb = b.Beta();
        double Bm = mu.Beta();

        const double Bm2 = Bm * Bm;
        Sx = (x0 * Bm - mu.P() * (1.0 - Bm * Bm)) / Bm2;
        Sy = (x0p / Bb - c * Sx) / s;

        w  = (Bm / Bb - c) / s;
        w_ = (-Bm / Bb - c) / s;

        Om2 = w * w + 1.0 - Bm * Bm;
        eps2 = (mW2 - mN2) * (1.0 - Bm * Bm);

        x1 = Sx - (Sx + w * Sy) / Om2;
        y1 = Sy - (Sx + w * Sy) * w / Om2;

        double Z2 = x1 * x1 * Om2 - (Sy - w * Sx) * (Sy - w * Sx) - (mW2 - x0 * x0 - eps2);

        // Calculate Z
        Z = std::sqrt(std::max(0.0, Z2));
    }

    // Extended rotation from F' to F coord.
    TMatrixD getK() const {
        TMatrixD K(4,4);
        K(0,0) = c; K(0,1) = -s; K(0,2) = 0; K(0,3) = 0;
        K(1,0) = s; K(1,1) =  c; K(1,2) = 0; K(1,3) = 0;
        K(2,0) = 0; K(2,1) =  0; K(2,2) = 1; K(2,3) = 0;
        K(3,0) = 0; K(3,1) =  0; K(3,2) = 0; K(3,3) = 1;
        return K;
    }

    // F coord. constraint on W momentum: ellipsoid.
    TMatrixD getA_mu_mat() const {
        TMatrixD A(4,4);
        const double B2 = mu.Beta() * mu.Beta();
        const double SxB2 = Sx * B2;
        double F = mW2 - x0 * x0 - eps2;

        A(0,0) = 1 - B2; A(0,1) = 0;    A(0,2) = 0; A(0,3) = SxB2;
        A(1,0) = 0;      A(1,1) = 1;    A(1,2) = 0; A(1,3) = 0;
        A(2,0) = 0;      A(2,1) = 0;    A(2,2) = 1; A(2,3) = 0;
        A(3,0) = SxB2;   A(3,1) = 0;    A(3,2) = 0; A(3,3) = F;

        return A;
    }

    // F coord. constraint on W momentum: ellipsoid
    TMatrixD A_b() const {
        const TMatrixD K = getK();
        const double B = b.Beta();

        TMatrixD A_b_(4, 4);
        A_b_(0,0) = 1 - B*B;  A_b_(0,1) = 0;      A_b_(0,2) = 0;      A_b_(0,3) = B*x0p;
        A_b_(1,0) = 0;        A_b_(1,1) = 1;      A_b_(1,2) = 0;      A_b_(1,3) = 0;
        A_b_(2,0) = 0;        A_b_(2,1) = 0;      A_b_(2,2) = 1;      A_b_(2,3) = 0;
        A_b_(3,0) = B*x0p;    A_b_(3,1) = 0;      A_b_(3,2) = 0;      A_b_(3,3) = mW2 - x0p*x0p;

        TMatrixD result = K * A_b_;
        TMatrixD KT(TMatrixD::kTransposed, K);
        result *= KT;
        return result;
    }

    // Rotation from F coord. to laboratory coord.
    TMatrixD getR_T() const {
        auto apply = [](const TMatrixD& R, const TVector3& vec) {
            TVector3 out;
            for (int i = 0; i < 3; ++i) {
                out[i] = R(i,0) * vec.X() + R(i,1) * vec.Y() + R(i,2) * vec.Z();
            }
            return out;
        };

        const TMatrixD Rz = Rotation(2, -mu.Phi());
        const TMatrixD Ry = Rotation(1, 0.5 * M_PI - mu.Theta());
        TVector3 rotated = apply(Ry, apply(Rz, b.Vect()));
        const TMatrixD Rx = Rotation(0, -std::atan2(rotated.Z(), rotated.Y()));

        TMatrixD RzT(TMatrixD::kTransposed, Rz);
        TMatrixD RyT(TMatrixD::kTransposed, Ry);
        TMatrixD RxT(TMatrixD::kTransposed, Rx);

        TMatrixD tmp = RyT * RxT;
        return RzT * tmp;
    }

    // H_tilde transformation from t=[c,s,1] to p_nu: F coord.
    TMatrixD getH_tilde() const {
        double _x1 = x1;
        double _y1 = y1;
        double p = mu.P();
        double _Z = Z;
        double _w = w;
        double _Om = std::sqrt(Om2);
        TMatrixD H_tilde(3, 3);
        H_tilde(0,0) =  _Z/_Om; H_tilde(0,1) = 0; H_tilde(0,2) = _x1-p;
        H_tilde(1,0) = _w*_Z/_Om; H_tilde(1,1) = 0; H_tilde(1,2) = _y1;
        H_tilde(2,0) = 0; H_tilde(2,1) = _Z; H_tilde(2,2) = 0;
        return H_tilde;
    }

    // Transformation of t=[c,s,1] to p_nu: lab coord.
    TMatrixD getH() const {
        TMatrixD result = getR_T() * getH_tilde();
        return result;
    }

    // Transformation of t=[c,s,1] to pT_nu: lab coord.
    TMatrixD getH_perp() const {
        TMatrixD h = getH();
        TMatrixD h_perp(3,3);
        h_perp.Zero();
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 3; ++j) {
                h_perp(i,j) = h(i,j);
            }
        }
        h_perp(2,2) = 1.0;
        return h_perp;
    }

    // Solution ellipse of pT_nu: lab coord.
    TMatrixD getN() const {
        // Invert H_perp safely, falling back to identity if singular
        TMatrixD Hp = getH_perp();
        TMatrixD HpInv(Hp);
        usedMatrixFallback_ = false;

        const double det = det3x3(Hp);
        if (!std::isfinite(det) || std::abs(det) < 1e-12) {
            usedMatrixFallback_ = true;
            TMatrixD fallback(3,3);
            fallback.UnitMatrix();
            return fallback;
        }

        HpInv.Invert();
        TMatrixD HpInvT(TMatrixD::kTransposed, HpInv);
        TMatrixD result = HpInvT * UnitCircle() * HpInv;
        return result;
    }

    bool usedMatrixFallback() const { return usedMatrixFallback_; }
};

// ---------- Class doubleNeutrinoSolution: finds the best double-neutrino solution for given b-jet and lepton momenta. ----------
class doubleNeutrinoSolution {
public:
    struct NuPair {
        std::array<double, 2> first;
        std::array<double, 2> second;
    };

    doubleNeutrinoSolution()
        : H1(3, 3),
          H2(3, 3),
          N1_(3, 3),
          N2_constraint_(3, 3),
          N2_nubar_(3, 3),
          usedMinimizerFallback_(false)
    {
        H1.Zero();
        H2.Zero();
        N1_.Zero();
        N2_constraint_.Zero();
        N2_nubar_.Zero();
    }

    doubleNeutrinoSolution(
        const TLorentzVector& b1,
        const TLorentzVector& b2,
        const TLorentzVector& l1,
        const TLorentzVector& l2,
        double met_x,
        double met_y)
        : H1(3, 3),
          H2(3, 3),
          N1_(3, 3),
          N2_constraint_(3, 3),
          N2_nubar_(3, 3),
          usedMinimizerFallback_(false)
    {
        H1.Zero();
        H2.Zero();
        N1_.Zero();
        N2_constraint_.Zero();
        N2_nubar_.Zero();
        const double mW = 80.385; // W mass in GeV
        const double mT = 172.5;  // Top mass in GeV

        auto try_pairing = [&](const TLorentzVector& B1, const TLorentzVector& B2,
                               const TLorentzVector& L1, const TLorentzVector& L2) {
            PairingResult result;

            nuana::nuSolutionSet ss1(B1, L1, mW, mT);
            nuana::nuSolutionSet ss2(B2, L2, mW, mT);

            TMatrixD H1tmp = ss1.getH();
            TMatrixD H2tmp = ss2.getH();
            result.H1.ResizeTo(H1tmp.GetNrows(), H1tmp.GetNcols());
            result.H2.ResizeTo(H2tmp.GetNrows(), H2tmp.GetNcols());
            result.H1 = H1tmp;
            result.H2 = H2tmp;

            TMatrixD V0(3, 3);
            V0.Zero();
            V0(0, 2) = met_x;
            V0(1, 2) = met_y;
            V0(2, 2) = 0.0;

            TMatrixD S = V0 - UnitCircle();
            TMatrixD ST(TMatrixD::kTransposed, S);

            TMatrixD N1 = ss1.getN();
            TMatrixD N2 = ss2.getN();

            TMatrixD n2 = ST * N2 * S;

            result.N1.ResizeTo(N1.GetNrows(), N1.GetNcols());
            result.N2_constraint.ResizeTo(n2.GetNrows(), n2.GetNcols());
            result.N2_nubar.ResizeTo(N2.GetNrows(), N2.GetNcols());
            result.N1 = N1;
            result.N2_constraint = n2;
            result.N2_nubar = N2;

            std::vector<TVectorD> intersections =
                nuana::intersections_ellipses(N1, n2).first;

            for (const auto& sol : intersections) {
                if (!satisfiesConic(sol, N1) || !satisfiesConic(sol, n2)) {
                    continue;
                }

                TVectorD nu1 = sol;
                TVectorD nu2 = S * sol;

                if (!satisfiesConic(nu2, N2)) {
                    continue;
                }

                double w1 = (nu1.GetNrows() > 2) ? nu1(2) : 1.0;
                double w2 = (nu2.GetNrows() > 2) ? nu2(2) : 1.0;

                if (!std::isfinite(w1) || !std::isfinite(w2) || std::abs(w1) <= 1e-12 || std::abs(w2) <= 1e-12) {
                    continue;
                }

                NuPair pair;
                pair.first = {nu1(0) / w1, nu1(1) / w1};
                pair.second = {nu2(0) / w2, nu2(1) / w2};

                if (isFinitePair(pair)) {
                    result.solutions.push_back(pair);
                }
            }

            if (result.solutions.empty() && kEnableMinimizerFallback) {
                TMatrixD es1 = ss1.getH_perp();
                TMatrixD es2 = ss2.getH_perp();

                TVectorD met_vec(3);
                met_vec(0) = met_x;
                met_vec(1) = met_y;
                met_vec(2) = 1.0;

                auto nus = [&](const TVectorD& ts) {
                    std::vector<TVectorD> momenta;
                    momenta.reserve(2);
                    for (int i = 0; i < 2; ++i) {
                        TVectorD vec(3);
                        vec(0) = std::cos(ts(i));
                        vec(1) = std::sin(ts(i));
                        vec(2) = 1.0;

                        TVectorD nu = (i == 0) ? es1 * vec : es2 * vec;
                        momenta.push_back(nu);
                    }
                    return momenta;
                };

                auto residuals = [&](const TVectorD& params) {
                    auto nu_vecs = nus(params);
                    TVectorD total = nu_vecs[0] + nu_vecs[1] - met_vec;
                    TVectorD res(2);
                    res(0) = total(0);
                    res(1) = total(1);
                    return res;
                };

                class ResidualsFunction : public ROOT::Math::IMultiGenFunction {
                public:
                    explicit ResidualsFunction(
                        const std::function<TVectorD(const TVectorD&)>& f)
                        : func(f) {}

                    unsigned int NDim() const override { return 2; }

                    double DoEval(const double* x) const override {
                        TVectorD params(2);
                        params(0) = x[0];
                        params(1) = x[1];
                        TVectorD res = func(params);
                        return res(0) * res(0) + res(1) * res(1);
                    }

                    ROOT::Math::IMultiGenFunction* Clone() const override {
                        return new ResidualsFunction(func);
                    }

                private:
                    std::function<TVectorD(const TVectorD&)> func;
                };

                std::unique_ptr<ROOT::Math::Minimizer> min(
                    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

                if (min) {
                    min->SetTolerance(1e-10);
                    min->SetPrecision(1e-12);
                    min->SetVariableStepSize(0, 0.01);
                    min->SetVariableStepSize(1, 0.01);

                    ResidualsFunction residualsFunc(residuals);
                    min->SetFunction(residualsFunc);

                    min->SetVariable(0, "t1", 0.0, 0.1);
                    min->SetVariable(1, "t2", 0.0, 0.1);

                    if (min->Minimize()) {
                        TVectorD ts(2);
                        ts(0) = min->X()[0];
                        ts(1) = min->X()[1];

                        auto fallbackSolutions = nus(ts);

                        NuPair pair;
                        pair.first = {fallbackSolutions[0](0), fallbackSolutions[0](1)};
                        pair.second = {fallbackSolutions[1](0), fallbackSolutions[1](1)};
                        result.solutions.push_back(pair);
                        result.usedMinimizerFallback = true;
                    }
                }
            }

            if (!result.solutions.empty()) {
                auto pairResidualSq = [&](const NuPair& pair) {
                    if (!std::isfinite(pair.first[0]) || !std::isfinite(pair.first[1]) ||
                        !std::isfinite(pair.second[0]) || !std::isfinite(pair.second[1])) {
                        return std::numeric_limits<double>::infinity();
                    }
                    double dx = (pair.first[0] + pair.second[0]) - met_x;
                    double dy = (pair.first[1] + pair.second[1]) - met_y;
                    double res = dx * dx + dy * dy;
                    return std::isfinite(res) ? res : std::numeric_limits<double>::infinity();
                };

                std::stable_sort(
                    result.solutions.begin(),
                    result.solutions.end(),
                    [&](const NuPair& a, const NuPair& b) {
                        return pairResidualSq(a) < pairResidualSq(b);
                    }
                );

            }
            return result;
        };

        PairingResult pairing1 = try_pairing(b1, b2, l1, l2);
        PairingResult pairing2 = try_pairing(b1, b2, l2, l1);

        double residual1 = metResidual(pairing1, met_x, met_y);
        double residual2 = metResidual(pairing2, met_x, met_y);

        if (residual1 <= residual2) {
            nunu_s = pairing1.solutions;
            H1.ResizeTo(pairing1.H1.GetNrows(), pairing1.H1.GetNcols());
            H2.ResizeTo(pairing1.H2.GetNrows(), pairing1.H2.GetNcols());
            H1 = pairing1.H1;
            H2 = pairing1.H2;
            N1_.ResizeTo(pairing1.N1.GetNrows(), pairing1.N1.GetNcols());
            N2_constraint_.ResizeTo(pairing1.N2_constraint.GetNrows(), pairing1.N2_constraint.GetNcols());
            N2_nubar_.ResizeTo(pairing1.N2_nubar.GetNrows(), pairing1.N2_nubar.GetNcols());
            N1_ = pairing1.N1;
            N2_constraint_ = pairing1.N2_constraint;
            N2_nubar_ = pairing1.N2_nubar;
            usedMinimizerFallback_ = pairing1.usedMinimizerFallback;
        } else {
            nunu_s = pairing2.solutions;
            H1.ResizeTo(pairing2.H1.GetNrows(), pairing2.H1.GetNcols());
            H2.ResizeTo(pairing2.H2.GetNrows(), pairing2.H2.GetNcols());
            H1 = pairing2.H1;
            H2 = pairing2.H2;
            N1_.ResizeTo(pairing2.N1.GetNrows(), pairing2.N1.GetNcols());
            N2_constraint_.ResizeTo(pairing2.N2_constraint.GetNrows(), pairing2.N2_constraint.GetNcols());
            N2_nubar_.ResizeTo(pairing2.N2_nubar.GetNrows(), pairing2.N2_nubar.GetNcols());
            N1_ = pairing2.N1;
            N2_constraint_ = pairing2.N2_constraint;
            N2_nubar_ = pairing2.N2_nubar;
            usedMinimizerFallback_ = pairing2.usedMinimizerFallback;
        }

        if (nunu_s.empty()) {
            NuPair empty_pair{{std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN()},
                              {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::quiet_NaN()}};
            nunu_s.push_back(empty_pair);
        }
    }

    std::vector<NuPair> get_nunu_s() const {
        return nunu_s;
    }

    const TMatrixD& getH1() const { return H1; }
    const TMatrixD& getH2() const { return H2; }
    const TMatrixD& getN1() const { return N1_; }
    const TMatrixD& getN2() const { return N2_constraint_; }
    const TMatrixD& getN2Constraint() const { return N2_constraint_; }
    const TMatrixD& getN2NuBar() const { return N2_nubar_; }

    size_t numSolutions() const { return nunu_s.size(); }

    bool isValid(size_t idx = 0) const {
        return idx < nunu_s.size() && isFinitePair(nunu_s[idx]);
    }

    bool hasValidSolution() const { return isValid(); }

    double nu1_px(size_t idx = 0) const {
        return idx < nunu_s.size() ? nunu_s[idx].first[0]
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double nu1_py(size_t idx = 0) const {
        return idx < nunu_s.size() ? nunu_s[idx].first[1]
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double nu2_px(size_t idx = 0) const {
        return idx < nunu_s.size() ? nunu_s[idx].second[0]
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    double nu2_py(size_t idx = 0) const {
        return idx < nunu_s.size() ? nunu_s[idx].second[1]
                                   : std::numeric_limits<double>::quiet_NaN();
    }

    bool usedMinimizerFallback() const { return usedMinimizerFallback_; }

    std::vector<double> allSolutionsFlat() const {
        std::vector<double> flat;
        flat.reserve(nunu_s.size() * 4);
        for (const auto& pair : nunu_s) {
            flat.push_back(pair.first[0]);
            flat.push_back(pair.first[1]);
            flat.push_back(pair.second[0]);
            flat.push_back(pair.second[1]);
        }
        return flat;
    }

private:
    TMatrixD H1, H2; // store the ellipse matrices of the selected pairing
    TMatrixD N1_;
    TMatrixD N2_constraint_;
    TMatrixD N2_nubar_;
    bool usedMinimizerFallback_ = false;

    struct PairingResult {
        std::vector<NuPair> solutions;
        TMatrixD H1;
        TMatrixD H2;
        TMatrixD N1;
        TMatrixD N2_constraint;
        TMatrixD N2_nubar;
        bool usedMinimizerFallback;

        PairingResult()
            : solutions(), H1(3, 3), H2(3, 3), N1(3, 3), N2_constraint(3, 3), N2_nubar(3, 3),
              usedMinimizerFallback(false)
        {
            H1.Zero();
            H2.Zero();
            N1.Zero();
            N2_constraint.Zero();
            N2_nubar.Zero();
        }
    };

    static bool isFinitePair(const NuPair& pair) {
        return std::isfinite(pair.first[0]) &&
               std::isfinite(pair.first[1]) &&
               std::isfinite(pair.second[0]) &&
               std::isfinite(pair.second[1]);
    }

    static double metResidual(const PairingResult& res, double met_x, double met_y) {
        if (res.solutions.empty()) {
            return std::numeric_limits<double>::infinity();
        }
        const auto& best = res.solutions.front();
        double sumx = best.first[0] + best.second[0];
        double sumy = best.first[1] + best.second[1];
        double residual = std::hypot(sumx - met_x, sumy - met_y);
        return residual;
    }

    std::vector<NuPair> nunu_s;
};

struct EventKinematics {
    TLorentzVector top1;
    TLorentzVector top2;
    double chel;
    double dphi_ttbar;
    double pdark;
    bool valid;

    EventKinematics()
        : top1(), top2(),
          chel(std::numeric_limits<double>::quiet_NaN()),
          dphi_ttbar(std::numeric_limits<double>::quiet_NaN()),
          pdark(std::numeric_limits<double>::quiet_NaN()),
          valid(false)
    {
    }
};

inline TLorentzVector makeMasslessNeutrino(double px, double py) {
    TLorentzVector nu;
    double energy = std::sqrt(px * px + py * py);
    nu.SetPxPyPzE(px, py, 0.0, energy);
    return nu;
}

inline TVector3 leptonDirectionInTopRest(const TLorentzVector& top, const TLorentzVector& lepton) {
    if (top.E() <= 0.0) {
        return TVector3();
    }

    TLorentzVector lep = lepton;
    lep.Boost(-top.BoostVector());
    TVector3 dir = lep.Vect();
    if (dir.Mag2() == 0.0) {
        return TVector3();
    }
    dir *= 1.0 / dir.Mag();
    return dir;
}

inline EventKinematics computeEventKinematics(const TLorentzVector& b1,
                                              const TLorentzVector& b2,
                                              const TLorentzVector& l1,
                                              const TLorentzVector& l2,
                                              double met_x,
                                              double met_y,
                                              const doubleNeutrinoSolution& solver,
                                              size_t idx = 0) {
    EventKinematics kin;

    if (!solver.isValid(idx)) {
        return kin;
    }

    TLorentzVector nu1 = makeMasslessNeutrino(solver.nu1_px(idx), solver.nu1_py(idx));
    TLorentzVector nu2 = makeMasslessNeutrino(solver.nu2_px(idx), solver.nu2_py(idx));

    kin.top1 = b1 + l1 + nu1;
    kin.top2 = b2 + l2 + nu2;

    TVector3 l1_rf = leptonDirectionInTopRest(kin.top1, l1);
    TVector3 l2_rf = leptonDirectionInTopRest(kin.top2, l2);
    kin.chel = l1_rf.Dot(l2_rf);

    kin.dphi_ttbar = std::abs(TVector2::Phi_mpi_pi(kin.top1.Phi() - kin.top2.Phi()));

    double residual_x = met_x - nu1.Px() - nu2.Px();
    double residual_y = met_y - nu1.Py() - nu2.Py();
    kin.pdark = residual_x * residual_x + residual_y * residual_y;

    kin.valid = true;
    return kin;
}
} // namespace nuana

// ---------- Aliases ----------
using nuana::nuSolutionSet;
using nuana::doubleNeutrinoSolution;
using nuana::makeMetCov;

#endif // MKSHAPESRDF_PROCESSOR_MODULES_DOUBLENEUTRINOSOLUTIONMACRO_H
