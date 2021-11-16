#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "Garfield/ComponentElmer.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char *argv[])
{
    TApplication app("app", &argc, argv);

    // tripple gem dimension (cm)
    const double pitch = 0.0140;
    const double tD = 0.0050; // Dielectric layer thickness
    
    const double dZp = 0.2; // Distance between pad plane(z = 0) and GEM3 layer
    const double dZ23 = 0.2; // betwen GEM 2 and GEM 3
    const double dZ12 = 0.2; // betwen GEM 1 and GEM 2
    const double dZu = 0.2; // betwen GEM1 and upper boundary plane
    
    const double gemX = 0.6;
    const double gemY = 0.6;
    // GEM voltages (V)
    const double dvg1 = 300;
    const double dvg2 = 300;
    const double dvg3 = 300;

    // electric field density (V/cm)
    const double eTrans = 1000;
    const double eDrift = 1000;
    const double eInduction = 1000;

    ComponentElmer fm(
        "tgemcell/mesh.header", "tgemcell/mesh.elements", "tgemcell/mesh.nodes",
        "tgemcell/dielectrics.dat", "tgemcell/tgemcell.result", "cm");
    fm.EnablePeriodicityX();
    fm.EnablePeriodicityY();
    fm.EnableConvergenceWarnings(false); // turn off convergence warning messages.
    MediumMagboltz gas;
    gas.SetComposition("Ar", 90, "iC4H10", 10);
    gas.SetPressure(760);
    gas.SetTemperature(293.15);
    gas.Initialise(true);
    gas.EnablePenningTransfer(0.51, 0., "Ar");
    const std::string path = std::getenv("GARFIELD_INSTALL");
    gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");
    
    const unsigned int nMaterials = fm.GetNumberOfMaterials();
    for (unsigned int i = 0; i < nMaterials; ++i)
    {
        const double eps = fm.GetPermittivity(i);
        if (eps == 1.)
            fm.SetMedium(i, &gas);
    }

    // plot equipotential lines of GEM1
    TCanvas *cf1 = new TCanvas("cf1", "Potential plot of GEM1", 600, 600); cf1->SetLeftMargin(0.16);
    TCanvas *cf2 = new TCanvas("cf2", "Potential plot of GEM2", 600, 600); cf2->SetLeftMargin(0.16);
    TCanvas *cf3 = new TCanvas("cf3", "Potential plot of GEM3", 600, 600); cf3->SetLeftMargin(0.16);
    TCanvas *cf4 = new TCanvas("cf4", "Potential plot", 600, 600); cf4->SetLeftMargin(0.16);

    ViewField vf;
    vf.SetComponent(&fm);
    vf.SetNumberOfContours(50);
    vf.SetPlaneXZ();

    vf.SetCanvas(cf1);
    vf.SetVoltageRange(
        -dZp * eInduction - dvg1 - dZ12 * eTrans - dvg2 - dZ23 * eTrans - 1.5 * dvg3,
        -dZp * eInduction - dvg1 - dZ12 * eTrans - dvg2 - dZ23 * eTrans + 1.5 * dvg3);
    vf.SetArea(-1.5 * pitch, dZp + dZ12 + dZ23 - 3 * tD, 1.5 * pitch, dZp + dZ12 + dZ23 + 3 * tD);
    vf.PlotContour();

    vf.SetCanvas(cf2);
    vf.SetVoltageRange(
        -dZp * eInduction - dvg1 - dZ12 * eTrans - 1.5 * dvg2,
        -dZp * eInduction - dvg1 - dZ12 * eTrans + 1.5 * dvg2);
    vf.SetArea(-1.5 * pitch, dZp + dZ12 - 3 * tD, 1.5 * pitch, dZp + dZ12 + 3 * tD);
    vf.PlotContour();

    vf.SetCanvas(cf3);
    vf.SetVoltageRange(
        -dZp * eInduction - 1.5 * dvg1,
        -dZp * eInduction + 1.5 * dvg1);
    vf.SetArea(-1.5 * pitch, dZp - 3 * tD, 1.5 * pitch, dZp + 3 * tD);
    vf.PlotContour();

    vf.SetCanvas(cf4);
    vf.SetVoltageRange(-dZp * eInduction - dvg1 - dZ12 * eTrans - dvg2 - dZ23 * eTrans - dvg3 - dZu*eDrift, 0.);
    vf.SetArea(-gemX/2, 0, gemX/2, dZp + dZ12 + dZ23 + dZu);
    vf.PlotContour();
    

    // simulating avalanche of an electron.
    ViewDrift driftView;

    Sensor sensor;
    sensor.AddComponent(&fm);
    sensor.SetArea(-gemX / 2, -gemY / 2, 0., gemX / 2, gemY / 2, dZp + dZ12 + dZ23 + dZu);

    AvalancheMicroscopic aval;
    aval.EnablePlotting(&driftView);
    aval.SetSensor(&sensor);
    // set the number of steps to be skipped when plotting.
    aval.SetCollisionSteps(100);

    // avalanching an electron
    const double x0 = 0.;
    const double y0 = 0.;
    const double z0 = dZp + dZ12 + dZ23 + dZu * 0.99; 
    const double t0 = 0.;
    const double e0 = 0.1;
    aval.AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    const unsigned int np = aval.GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (unsigned int j = 0; j < np; ++j) {
      aval.GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                  xe2, ye2, ze2, te2, e2, status);
      // drift.DriftIon(xe1, ye1, ze1, te1);
      // drift.GetIonEndpoint(0, xi1, yi1, zi1, ti1, xi2, yi2, zi2, ti2, status);
    }
    std::cout << "The size of avalanche : " << ne << "\n";

    TCanvas *cd1 = new TCanvas("cd1", "Drift line of GEM1", 600, 600);cd1->SetLeftMargin(0.16);
    TCanvas *cd2 = new TCanvas("cd2", "Drift line of GEM2", 600, 600);cd2->SetLeftMargin(0.16);
    TCanvas *cd3 = new TCanvas("cd3", "Drift line of GEM3", 600, 600);cd3->SetLeftMargin(0.16);
    TCanvas *cd4 = new TCanvas("cd4", "Drift line", 600, 600);cd4->SetLeftMargin(0.16);

    ViewFEMesh meshView;
    meshView.SetComponent(&fm);
    // x-z projection.
    meshView.SetPlaneXZ();
    // Set the color of the kapton.
    meshView.SetFillMesh(true);
    meshView.SetColor(0, kGray);
    meshView.SetColor(3, kYellow + 3);
    meshView.SetColor(6, kYellow + 3);
    meshView.SetColor(9, kYellow + 3);
    meshView.SetViewDrift(&driftView);
    meshView.EnableAxes();

    meshView.SetCanvas(cd1);
    meshView.SetArea(-1.5 * pitch, dZp + dZ12 + dZ23 - 3 * tD, 1.5 * pitch, dZp + dZ12 + dZ23 + 3 * tD);
    meshView.Plot();

    meshView.SetCanvas(cd2);
    meshView.SetArea(-1.5 * pitch, dZp + dZ12 - 3 * tD, 1.5 * pitch, dZp + dZ12 + 3 * tD);
    meshView.Plot();

    meshView.SetCanvas(cd3);
    meshView.SetArea(-1.5 * pitch, dZp - 3 * tD, 1.5 * pitch, dZp + 3 * tD);
    meshView.Plot();

    meshView.SetCanvas(cd4);
    meshView.SetArea(-gemX/2, 0, gemX/2, dZp + dZ12 + dZ23 + dZu);
    meshView.Plot();

    app.Run(true);
}
