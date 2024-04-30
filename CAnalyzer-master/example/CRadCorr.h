#ifndef C_RAD_CORR_H
#define C_RAD_CORR_H

#include <cmath>
#include <vector>
#include <string>
#include "ConfigObject.h"
#include "CExpData.h"
#include "CModelWrapper.h"


class CRadCorr : public ConfigObject
{
public:
    CRadCorr(const std::string &path = "");
    virtual ~CRadCorr();

    void Configure(const std::string &path);
    void RadiativeCorrection(CExpData &exp_data, int iters = 0);
    void Radiate(CExpData &exp_data);
    void Initialize(const CExpData &exp_data, bool radiate = false);
    bool SanityCheck(const CExpData &exp_data);
    const CModelWrapper &GetModel() {return model;}

private:
    void radcor(CExpData::DataSet &dset, bool verbose = true);
    void xyrad2d(CExpData::DataSet &dset, bool verbose = true);
    double fes(double Es);
    double fep(double Ep);
    double int_es(double Es);
    double int_esep(double Ep, double Es);
    double int_ep(double Ep);
    double int_esdp(double Es);
    double get_cxsn(double E0, double Eb);
    void init_model(const CExpData::DataSet &ref_set);
    void scale_model(const CExpData &exp_data, bool born_level);
    template<typename T>
    void number_operation(const std::string &key, T &val);

    // some inlines
    inline CExpData::DataSet get_ref_set(const CExpData &exp_data, bool model);
    inline void spectrum_init(const CExpData::DataSet &dset);
    inline void point_init(const CExpData::DataPoint &point);
    inline double __Ep_max(double Es);
    inline double __Es_min(double Ep);
    inline double __phi(double x);
    inline double __eta(double Z);
    inline double __Q2(double E, double Epr);
    inline double __log_Q2m2(double E, double Epr);
    inline double __F_bar(double E, double Epr, double gamma_t);
    inline double __btr(double E, double Epr);
    inline double __I(double E0, double E, double xi, double bt);
    inline double __XI_Stein(double radl);

private:
    bool internal_RC, external_RC, user_defined_XI, peak_approx, use_model;
    int n_sim, n_sim_2d;
    double iter_prec;
    double target_Z, target_A, target_M;
    double theta, sin2, cos2;

    // parameters that will be shared between different functions
    double F_mott, Schwinger, delta, delta1, delta2, Bz; // for whole data sets
    double Es, BTB, BTA, XIB, XIA, GAMT;                 // for each spectrum
    double Ep, R, BTR, Epmin, Epmax, Esmin, Esmax;       // for each data point

    // model related
    CModelWrapper model;
};

#endif
