#ifndef TOPNESS_PRODUCER_CC
#define TOPNESS_PRODUCER_CC

#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

using namespace ROOT;
using namespace ROOT::VecOps;

const int DMAX = 4;

struct my_func {
    virtual double operator()(double x[], int d) = 0;
    virtual ~my_func() = default;
};

class my_simplex {

public:
    int d;
    double alpha, beta, gamma;
    my_func *f;

    double xstart[DMAX*(DMAX+1)];
    double x[DMAX*(DMAX+1)];

    double xh[DMAX];
    double xl[DMAX];

    double y[DMAX+1];
    double yl, ynh, yh;

    double xCentroid[DMAX];
    double xReflect[DMAX];
    double xExpand[DMAX];
    double xContract[DMAX];

    double yReflect, yExpand, yContract;

    int imin, imax, inmax;

    my_simplex(int dd,double aalpha,double bbeta,double ggamma,my_func *ff)
        : d(dd),alpha(aalpha),beta(bbeta),gamma(ggamma),f(ff){}

    void my_SetUp(double xin[]){
        int D=d*(d+1);

        std::copy(xin,xin+D,xstart);

        for(int i=0;i<D;i++)
            x[i]=xstart[i];

        for(int i=0;i<d+1;i++){
            double xi[DMAX];
            std::copy(x+d*i,x+d*i+d,xi);
            y[i]=(*f)(xi,d);
        }
    }

    void set_y(){
        for(int i=0;i<d+1;i++){
            double xi[DMAX];
            std::copy(x+d*i,x+d*i+d,xi);
            y[i]=(*f)(xi,d);
        }
    }

    void find_max(){
        for(int i=0;i<d+1;i++)
            if(y[i]>y[imax]) imax=i;

        if(imax==1) inmax=0;
        for(int i=0;i<d+1;i++){
            if(i==imax) continue;
            if(y[i]>y[inmax]) inmax=i;
        }

        yh=y[imax];
        ynh=y[inmax];

        for(int j=0;j<d;j++)
            xh[j]=x[d*imax+j];
    }

    void find_min(){
        imin=0;
        for(int i=0;i<d+1;i++)
            if(y[i]<y[imin]) imin=i;
        
	yl=y[imin];
        for(int j=0;j<d;j++)
            xl[j]=x[d*imin+j];
    }

    void my_Centroid(int h){
        for(int j=0;j<d;j++){
            xCentroid[j]=0;
            for(int i=0;i<d+1;i++){
                if(i==h) continue;
                xCentroid[j]+=x[d*i+j];
            }
        }
    }

    void my_Reflection(){
        for(int i=0;i<d;i++)
            xReflect[i]=xCentroid[i]*(1+alpha)/d-alpha*xh[i];
        yReflect=(*f)(xReflect,d);
    }

    void my_Expansion(){
        for(int j=0;j<d;j++)
            xExpand[j]=xCentroid[j]*(1-gamma)/d+gamma*xReflect[j];
        yExpand=(*f)(xExpand,d);
    }

    void my_Contraction(){
        for(int j=0;j<d;j++)
            xContract[j]=xCentroid[j]*(1-beta)/d+beta*xh[j];
        yContract=(*f)(xContract,d);
    }

    void replace_all(){
        for(int i=0;i<d+1;i++){
            if(i==imin) continue;
            for(int j=0;j<d;j++)
                x[d*i+j]=0.5*(x[d*i+j]+xl[j]);
        }
    }
};

class my_Nelder_Mead{

public:
    int d,Ntry,n;
    double eps,Deltastep;

    my_func *f;
    my_simplex simplex;

    bool convergeYes;
    double yfinal;
    double xfinal[DMAX*(DMAX+1)];

    my_Nelder_Mead(int dd,double alpha,double beta,double gamma,
                   int NNtry,double eeps,double DDeltastep,int nn,my_func *ff)
        : d(dd),Ntry(NNtry),eps(eeps),Deltastep(DDeltastep),n(nn),
          f(ff),simplex(dd,alpha,beta,gamma,f){}

    bool one_cycle(my_simplex *s){
        s->my_Centroid(s->imax);
        s->my_Reflection();
        if(s->yReflect<=s->yl){
            s->my_Expansion();
            if(s->yExpand<s->yl){
                std::copy(s->xExpand,s->xExpand+d,s->x+d*s->imax);
                s->set_y();
                return false;
            }

            else{
                std::copy(s->xReflect,s->xReflect+d,s->x+d*s->imax);
                s->set_y();
                return true;
            }
        }
        else if(s->yReflect>=s->ynh){
            if(s->yReflect<s->yh){
                std::copy(s->xReflect,s->xReflect+d,s->x+d*s->imax);
                s->set_y();
                s->find_max();
            }
            s->my_Contraction();

            if(s->yContract<s->yh){
                std::copy(s->xContract,s->xContract+d,s->x+d*s->imax);
                s->set_y();
                return false;
            }

            else{
                s->replace_all();
                s->set_y();
                return false;
            }
        }

        else{
            std::copy(s->xReflect,s->xReflect+d,s->x+d*s->imax);
            s->set_y();
            return true;
        }
    }

    bool find_global_min(double xin[]){
        yfinal=1e13;
        simplex.imax=0;
        simplex.inmax=1;
        simplex.imin=2;
        simplex.my_SetUp(xin);
        simplex.find_max();
        simplex.find_min();

        convergeYes=false;

        for(int i=0;i<Ntry;i++){
            one_cycle(&simplex);

            simplex.find_max();
            simplex.find_min();

            double ynewmin=simplex.yl;
            double ynewmax=simplex.yh;

            if(std::abs(ynewmax-ynewmin)/
               (std::abs(ynewmax)+std::abs(ynewmin)+eps)<eps){
                convergeYes=true;

                if(ynewmin<yfinal){

                    std::copy(simplex.x,
                              simplex.x+d*(d+1),
                              xfinal);
                    yfinal=ynewmin;
                    break;
                }
            }
        }

        return convergeYes;
    }
};

struct TopnessObjective: my_func{

    TLorentzVector b1,b2,lep;
    double nu_px,nu_py;

    TopnessObjective(const TLorentzVector &ib1,
                     const TLorentzVector &ib2,
                     const TLorentzVector &ilep,
                     double inu_px,
                     double inu_py)
        :b1(ib1),b2(ib2),lep(ilep),nu_px(inu_px),nu_py(inu_py){}

    double operator()(double x[],int d){

        double pWx=x[0];
        double pWy=x[1];
        double pWz=x[2];
        double nu_pz=x[3];

        const double mW=80.379;
        const double mt=172.5;

        const double aW=5.0;
        const double at=15.0;
        const double aCM=1000.0;

        TLorentzVector nu;
        nu.SetPxPyPzE(nu_px,nu_py,nu_pz,
        std::sqrt(nu_px*nu_px+nu_py*nu_py+nu_pz*nu_pz));

        TLorentzVector W;
        W.SetPxPyPzE(pWx,pWy,pWz,
        std::sqrt(pWx*pWx+pWy*pWy+pWz*pWz+mW*mW));

        TLorentzVector top1=b1+lep+nu;
        TLorentzVector top2=b2+W;
        TLorentzVector sum=b1+b2+lep+nu+W;

        double S=
        std::pow(mW*mW-W.M2(),2)/std::pow(aW,4)
       +std::pow(mt*mt-top1.M2(),2)/std::pow(at,4)
       +std::pow(mt*mt-top2.M2(),2)/std::pow(at,4)
       +std::pow(4*mt*mt-sum.M2(),2)/std::pow(aCM,4);

        return S;
    }
};

double topness_producer(

int nCleanJet,
RVecF CleanJet_pt,
RVecF CleanJet_eta,
RVecF CleanJet_phi,
RVecF CleanJet_mass,

int nLep,
RVecF Lep_pt,
RVecF Lep_eta,
RVecF Lep_phi,
RVecI Lep_pdgId,

float PuppiMET_pt,
float PuppiMET_phi,

RVecF Jet_btagger,
float bAlgo_WP
){
    if(nLep<1||nCleanJet<2) return NAN;

    TLorentzVector lep;

    double lepMass=(std::abs(Lep_pdgId[0])==13)?
                    0.105658:0.000511;

    lep.SetPtEtaPhiM(
        Lep_pt[0],
        Lep_eta[0],
        Lep_phi[0],
        lepMass);

    std::vector<int>bjets;

    for(size_t i=0;i<CleanJet_pt.size();i++){
        if(CleanJet_pt[i]>30 &&
           std::abs(CleanJet_eta[i])<2.5 &&
           Jet_btagger[i]>bAlgo_WP)
            bjets.push_back(i);
    }

    std::vector<std::pair<int,int>> jetPairs;

    if (bjets.size() >= 2) {
        for(size_t i=0;i<bjets.size();i++){
            for(size_t j=i+1;j<bjets.size();j++){
                jetPairs.push_back({bjets[i], bjets[j]});
            }
        }

    }
    else if (bjets.size() == 1) {
        std::vector<int> untagged;
        for (size_t i = 0; i < CleanJet_pt.size(); ++i) {
            if ((int)i == bjets[0]) continue;
            if (CleanJet_pt[i] > 30.0)
                untagged.push_back(i);
        }
        std::sort(untagged.begin(), untagged.end(),
            [&](int a, int b) { return CleanJet_pt[a] > CleanJet_pt[b]; });
        for (int i = 0; i < std::min(2,(int)untagged.size()); ++i) {
            jetPairs.push_back({bjets[0], untagged[i]});
        }
    }
    else {
        return NAN;
    }

    double nu_px=PuppiMET_pt*cos(PuppiMET_phi);
    double nu_py=PuppiMET_pt*sin(PuppiMET_phi);

    double bestS = 1e13;
    double xfinal[DMAX];

    for(auto &pair : jetPairs){

        TLorentzVector b1,b2;

        int j1 = pair.first;
        int j2 = pair.second;

        b1.SetPtEtaPhiM(
            CleanJet_pt[j1],
            CleanJet_eta[j1],
            CleanJet_phi[j1],
            CleanJet_mass[j1]);

        b2.SetPtEtaPhiM(
            CleanJet_pt[j2],
            CleanJet_eta[j2],
            CleanJet_phi[j2],
            CleanJet_mass[j2]);

        double alpha=1.0;
        double beta=0.5;
        double gamma=2.0;

        double eps=0.000002;
        double Deltastep=20.;

        int Ntry=100000;
        int Nattempts=100;
        int n=3000;

        srand((unsigned)time(0));

        double ybest1 = 1e13;
        double ybest2 = 1e13;
        double xbest1[DMAX] = {1000.,1000.,1000.,1000.};
        double xbest2[DMAX] = {1000.,1000.,1000.,1000.};

        double edir[DMAX*DMAX] = {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        };

        for(int attempt = 0; attempt < Nattempts; attempt++) {

            double xin[DMAX*(DMAX+1)];
            double xstart[DMAX*(DMAX+1)];
 
            for(int i = 0; i < DMAX+1; i++) {
                for(int j = 0; j < DMAX; j++) {

                    double start = ((rand() % n) / (double)n - 0.5); // [-0.5,0.5]

                    if(i == 0) {
                        xin[i*DMAX + j] = 8000.0 * start;
                    } else {
                        xin[i*DMAX + j] = xin[j] + Deltastep * edir[DMAX*(i-1) + j];
                    }
                }
                // Copy current simplex point into xstart
                std::copy(xin + i*DMAX, xin + i*DMAX + DMAX, xstart + i*DMAX);
            }

            // First combination (b1, b2)
            TopnessObjective objective(b1,b2,lep,nu_px,nu_py);
            my_Nelder_Mead nm1(DMAX, alpha, beta, gamma, Ntry, eps, Deltastep, n, &objective);
            bool converge1 = nm1.find_global_min(xstart);

            if(converge1 && nm1.yfinal < ybest1) {
                ybest1 = nm1.yfinal;
                std::copy(nm1.xfinal, nm1.xfinal + DMAX, xbest1);
            }

            // Second combination (b2, b1)
            TopnessObjective objective2(b2, b1, lep, nu_px, nu_py);
            my_Nelder_Mead nm2(DMAX, alpha, beta, gamma, Ntry, eps, Deltastep, n, &objective2);
            bool converge2 = nm2.find_global_min(xstart);
            if(converge2 && nm2.yfinal < ybest2) {
                ybest2 = nm2.yfinal;
                std::copy(nm2.xfinal, nm2.xfinal + DMAX, xbest2);
            }
        }
        // Choose the absolute best
        if(ybest1 < ybest2) {
            bestS = ybest1;
            std::copy(xbest1, xbest1 + DMAX, xfinal);
        } else {
            bestS = ybest2;
            std::copy(xbest2, xbest2 + DMAX, xfinal);
        }
    }
    if (bestS >= 1e13) {
//    std::cout << "Minimization failed for all jet pairs!" << std::endl;
    return NAN;
    }

    return std::log(bestS);
}

#endif
