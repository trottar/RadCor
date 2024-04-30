#include "CExpData.h"
#include "CRadCorr.h"

using namespace std;

// void rad_corr(const string &data_conf = "configs/data_sets_9deg.conf")
//void rad_corr(const string &data_conf = "configs/data_sets_test.conf", const string &output = "output/radcor_out.dat")
void rad_corr(const string &data_conf, const string &output)
{
    CRadCorr rad_cor;
    rad_cor.Configure("configs/rad_corr.conf");

    cout << "Internal RC is ";
    if(rad_cor.GetConfig<bool>("Internal RC"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "External RC is ";
    if(rad_cor.GetConfig<bool>("External RC"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "User Defined XI is ";
    if(rad_cor.GetConfig<bool>("User Defined XI"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "Peaking Approximation is ";
    if(rad_cor.GetConfig<bool>("Energy Peaking Approximation"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    CExpData data;
    data.ReadConfigFile(data_conf);

//    rad_cor.Radiate(data);
    rad_cor.RadiativeCorrection(data);

    data.SaveResult(output);
}

// void radiate(const string &data_conf = "configs/dxs_model_6deg.conf")
void radiate(const string &data_conf = "configs/data_sets_test.conf", const string &output = "output/radiated_model.dat")
{
    CRadCorr rad_cor;
    rad_cor.Configure("configs/rad_corr.conf");

    cout << "Internal RC is ";
    if(rad_cor.GetConfig<bool>("Internal RC"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "External RC is ";
    if(rad_cor.GetConfig<bool>("External RC"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "User Defined XI is ";
    if(rad_cor.GetConfig<bool>("User Defined XI"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    cout << "Peaking Approximation is ";
    if(rad_cor.GetConfig<bool>("Energy Peaking Approximation"))
        cout << "ON" << endl;
    else
        cout << "OFF" << endl;

    CExpData data;
    data.ReadConfigFile(data_conf);

    rad_cor.Radiate(data);

    data.SaveResult(output);
}

void radiate_all() {
    string data_conf, output;
    vector<string> thetas = {"11", "18", "30"};
    vector<string> pols = {"unpol", "long", "trans"};
    for (auto theta : thetas) {
        for (auto pol : pols) {
            data_conf = "configs/data_sets_radiate_" + theta + "deg_" + pol + ".conf";
            output = "output/radiated_model_" + theta + "deg_" + pol + ".dat";
            cout << endl << "Radiating for theta = " + theta + " and Pol = " + pol << endl;
            cout << "Using data_config file: " + data_conf << endl;
            for (int i = 0; i < 20; i++) cout << "-";
            cout << endl;
            radiate(data_conf, output);
            for (int i = 0; i < 20; i++) cout << "-";
            cout << endl;
            cout << "Radiated model output file: " + output << endl;
        }
    }
}
