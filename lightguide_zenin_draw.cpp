//check pos dependence
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>

#include "/Users/az/Software/root/ROBAST/include/AOpticsManager.h"
#include "/Users/az/Software/root/ROBAST/include/AFocalSurface.h"
#include "/Users/az/Software/root/ROBAST/include/ALens.h"
#include "/Users/az/Software/root/ROBAST/include/AMirror.h"
#include "/Users/az/Software/root/ROBAST/include/AGeoBezierPgon.h"
#include "/Users/az/Software/root/ROBAST/include/AGeoBezierPcon.h"
#include "/Users/az/Software/root/ROBAST/include/AGeoWinstonConePoly.h"
#include "/Users/az/Software/root/ROBAST/include/AGeoWinstonCone2D.h"
#include "/Users/az/Software/root/ROBAST/include/ARayShooter.h"

#include "/Users/az/Software/root/include/TGeoBBox.h"
#include "/Users/az/Software/root/include/TGeoSphere.h"
#include "/Users/az/Software/root/include/TGeoMatrix.h"
#include "/Users/az/Software/root/include/TCanvas.h"
#include "/Users/az/Software/root/include/TFile.h"
#include "/Users/az/Software/root/include/TGeoCompositeShape.h"
#include "/Users/az/Software/root/include/TAxis.h"
#include "/Users/az/Software/root/include/TLegend.h"
#include "/Users/az/Software/root/include/TGraph.h"
#include "/Users/az/Software/root/include/TGeoPgon.h"
#include "/Users/az/Software/root/include/TGeoParaboloid.h"
#include "/Users/az/Software/root/include/TRandom.h"
#include "/Users/az/Software/root/include/TMultiGraph.h"
#include "/Users/az/Software/root/include/TH2D.h"
#include "/Users/az/Software/root/include/TGaxis.h"
#include "/Users/az/Software/root/include/TGeoManager.h"
#include "/Users/az/Software/root/include/TThread.h"
#include "/Users/az/Software/root/include/TGraphErrors.h"

#define PI 3.14159265
// define useful units
static const Double_t cm = AOpticsManager::cm();
static const Double_t mm = AOpticsManager::mm();
static const Double_t um = AOpticsManager::um();
static const Double_t nm = AOpticsManager::nm();
static const Double_t  m = AOpticsManager::m();


//control section
const int mode = 3; //0 - random stuff/checks/plots; 1 - collection eff; 3 - pde comp; 2 - obsolete; 4 - absorption comp; 5 SNR/NSB
int cone = 0; //for plotting
const int spectrum = 0; //for now should both have the same value
const int cam_weight = spectrum; //for now should both have the same value
const int nsb = 0;
const int filter = 0;
const string pde_type = "old";

const int usePosWeight = 1;
const int useAnodePosWeight = 0;
const int var_ref = 1; //variable reflectance
const int useAngularWeight = 1;
const int useOdistr = 0;

const int usemppc = 1;
unsigned int coating = 0;
unsigned int epoxy = 0;
unsigned int absorption = 1;
unsigned int hamacut = 1;
const Double_t inc_angle = 0;

const int draw = 0;
const int test = 1;
const int ph_number = 20000;
const int intlambda = 310;
Double_t lambda = intlambda*nm; //wavelength

const Double_t pmt_fold = 1.0;
const Double_t sipm_fold = 1.0;

const Double_t start = 0;
const Double_t finish = 40;
const int steps = 41;

//filter
const Double_t wl_left = 490;
const Double_t wl_right = 560;
const Double_t eff = 0.87; //efficiency at wl < cutoff






//define parameters of the cone
//const Double_t kHexSize = 49.8*mm;  //flat to flat
//const Double_t kSquareSize = 25.0*mm; //MPPC size
//const Double_t kSideLength = (kHexSize/2.)/TMath::Cos(30*PI/180);
//const Double_t kConeHeight = 6.8*cm;
//const Double_t kEdgeThickness = 1*mm;
//const Double_t kOuterHexSize = kHexSize+kEdgeThickness;
//const Double_t kOuterSideLength = (kOuterHexSize/2.)/TMath::Cos(30*PI/180);


//mppc param
const Double_t kGap = 0.2*mm;
const Double_t kActiveWidth = 3.0*mm;
const Double_t kChipThickness = 0.1*mm;
const Double_t kChipWidth = 3.1*mm;
const Double_t kWidth = 25.6*mm;
const Double_t kCellWidth = 50*um;
const Double_t kSilThickness = 1.5*nm; //absorption layer
const Double_t kCellThickness = 1*um;
const Double_t kResinThickness = 0.1*mm;
const Double_t SiOThickness = 1000*nm;
const Double_t SiNThickness = 1000*nm;
const Double_t WallThickness = 2.*um;
const Double_t CellScaleFactor = 0.86; //root of 74%

const Double_t kMppcThickness = kSilThickness + kChipThickness +kResinThickness + SiOThickness + SiNThickness;

const Double_t kMppc_qe = 0.998; //48% PDE @ 405nm
const Double_t kPmt_ce = 0.95; //collection efficiency

Double_t ResinThickness = 0;
//Double_t mppc_offset = 3.5 * mm;

//parameters
//Double_t stepH = 0.05, stepPara = 0.05; //steps for optimization cycles
//Double_t paraF_r = 0.67, paraF_z = 0.86; //flat section parameters
//Double_t paraS_r = 0.3, paraS_z = 0.6; //side section parameters
//Double_t paraC_r = 0.3, paraC_z = 0.6; //corner section parameters
//Double_t paraH = kConeHeight; //height parameters
Double_t phfrac = 0; //fraction of focused photons

//O cone parameters
Double_t fRin = (49.1/2)*mm;
Double_t fRout = (24.4/2)*mm;
Double_t fRpmt = 21*mm;
Double_t height = 0;

//Double_t fRin = 25*mm;
//Double_t fRout = 11*mm;
//Double_t fRpmt = 21*mm;
//Double_t height = 68*mm;
//Double_t height = 0;

const Double_t kSideLength = (fRin)/TMath::Cos(30*PI/180);
//const Double_t kConeHeight = 6.8*cm;
const Double_t kEdgeThickness = 1*mm;
const Double_t kOuterHexSize = fRin*2+2*kEdgeThickness;
const Double_t kOuterSideLength = (kOuterHexSize/2.)/TMath::Cos(30*PI/180);


//graph scale function
void rescaleX(TGraph* gra, Double_t scale)
{
    int N = gra->GetN();
    Double_t* X = gra->GetX();
    for(int i = 0; i < N; i++)
        X[i]*=scale;
    gra->GetHistogram()->Delete();
    gra->SetHistogram(0);
    return;
}

void rescaleY(TGraph* gra, Double_t scale)
{
    int N = gra->GetN();
    Double_t* Y = gra->GetY();
    for(int i = 0; i < N; i++)
        Y[i]*=scale;
    gra->GetHistogram()->Delete();
    gra->SetHistogram(0);
    return;
}

void fold(TGraphErrors* gra, double offset = 0.)
{
    int N = gra->GetN();
    Double_t* Y = gra->GetY();
    Double_t* X = gra->GetX();
    Double_t* eY = gra->GetEY();
    Double_t* eX = gra->GetEX();
    for(int i = 0; i < N; i++)
    {
        gra->SetPoint(i, fabs(X[i] - offset), Y[i]);
        gra->SetPointError(i, eX[i], eY[i]);
    }
    return;
}

double filter_gauss(double x)
{
    double rlim = wl_right * nm;
    double llim = wl_left * nm;
    double mean = (rlim + llim) / 2.;
    double std = (rlim - llim) / 3;
    TF1* erf_in = new TF1("erf_in", "1/(TMath::Sqrt(3.1415 * [1]**2))*TMath::Exp(-((x-[0])/[1])**2)", mean - 4 * std, mean + 4 * std);
    erf_in->SetParameters(mean, std);
    if(x <= mean - 4 * std)
        return eff;
    else if(x >= mean - 4 * std && x <= mean + 4 * std)
        return (1 - erf_in->Integral(-4 * std, x)) * eff;
    else if(x >= mean + 4 * std)
        return 0;
}

//graph functions
TGraph* QEGraph2()
{
  TGraph* graQE = new TGraph;

  TFile f("/Users/az/CTA/work/lightguide/170411LightGuideROBAST/can_angle.root");

  TGraph* graph = (TGraph*)((TCanvas*)f.Get("can_angle"))->GetPrimitive("graph_ave3");

  Int_t n = graph->GetN();
  Double_t ymax = TMath::MaxElement(n, graph->GetY());

  for(Int_t i = 0; i < n; i++){
    Double_t y = graph->GetY()[i];
    Double_t x = graph->GetX()[i];
    graQE->SetPoint(i, x*TMath::DegToRad(), y/ymax);
  } // i

  delete graph;

  return graQE;
}

TGraph* plotgraph()
{
  //TGraph* plogra = new TGraph;

  TFile f("/Users/az/CTA/work/lightguide/170411LightGuideROBAST/resin_abs.root");
//    TFile f("/Users/az/CTA/work/lightguide/170411LightGuideROBAST/point_465nm_0deg_v2_low.root");

  TGraph* graph = (TGraph*)f.Get("graph");
  rescaleY(graph, 100);

return graph;
}

TGraph* QEGraph1()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/qe_R11920-100-02.dat", "%lf\t%lf");
    rescaleX(g,nm);
    return g;
}

//TGraph* hamaPDEGraph()
//{
//    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/pde_S14521-8648.dat", " %lf\t%lf");
//    rescaleX(g,nm);
//    rescaleY(g,0.01);
//    return g;
//}

TGraph* hamaPDEGraph75SC()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/pde_S14521-8648.dat", " %lf\t%lf");
    rescaleX(g,nm);
    rescaleY(g,0.01*0.74/0.82*0.893*1.125); // rescale to percent unit, channel fill factor ratio, whole mppc fill factor, PDE (3v->7v OV)
    return g;
}

TGraph* hamaPDEGraph()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/pde_S14520-3050.dat", " %lf\t%lf");
    rescaleX(g,nm);
    rescaleY(g,0.01);
    return g;
}

TGraph* hamaPDEGraphSC()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/pde_S14520-3050.dat", " %lf\t%lf");
    rescaleX(g,nm);
    rescaleY(g,0.01*0.893*1.125); // rescale to percent unit, channel fill factor ratio, whole mppc fill factor
    return g;
}

TGraph* OldHamaPDEGraphSC()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/hamamatsu.dat", " %lf\t%lf");
    rescaleX(g,nm);
    rescaleY(g,0.01*0.893*1.13); // rescale to percent unit, whole mppc fill factor
    return g;
}

TGraph* hamaQEGraph()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/qe_R12992-100-05_new.dat", "%lf,%lf");
    rescaleX(g,nm);
    rescaleY(g, 0.01);
    return g;
}

TGraph* siAbslenGraph()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/si_abs.dat", "%lf\t%lf");
    rescaleX(g,nm);
    rescaleY(g,cm);
    return g;
}

TGraph* CherenkovGraph()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/cherenkov.dat", "%lf\t%lf");
    rescaleX(g,nm);
    return g;
}

TGraph* NSBGraph()
{
    TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/NSB.dat", "%lf\t%lf");
    rescaleX(g,nm);
//    rescaleY(g,0.1);
    return g;
}


TGraph* AbslenGraph()
{
    TGraph* g;
    TFile f("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/resin_abs.root");

    if(!epoxy)
    {
        g = (TGraph*)((TCanvas*)f.Get("can3"))->GetPrimitive("silicone");
    }
    else
    {
        g = (TGraph*)((TCanvas*)f.Get("can3"))->GetPrimitive("epoxy");
    }
    rescaleX(g,nm);
    rescaleY(g,mm);
    return g;
}

Double_t decayfun_new(Double_t x)
{
    if(hamacut == 0)
        return 1.;
    else
    {
        TGraph *g;
        if(pde_type == "50")
            g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/scale_hama50.txt", "%lf\t%lf");
        else if(pde_type == "75")
            g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/scale_hama75.txt", "%lf\t%lf");
        else if(pde_type == "old")
            g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/scale_hama_old.txt", "%lf\t%lf");
        else
            cout << "Error in pde selection" << endl;
        rescaleX(g,nm);
        return g->Eval(x, 0, "S");;
    }
}

Double_t decayfun(Double_t x)
{
    if(hamacut == 0)
        return 1.;
    else
    {
        TGraph *g = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/scale_hama_old.txt", "%lf\t%lf");
        rescaleX(g,nm);
        return g->Eval(x, 0, "S");;
    }
}

Double_t CamWeight(Double_t angle)
{
    if(angle < 0)
        angle*=(-1);
    if(angle > 40.)
        return 0;
    else
        return angle/40. * 10;
}

TH1D* CamWeightNSB()
{
    TFile* f = new TFile("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/nsb_distr_3.root");
    TCanvas* c = (TCanvas*) f->Get("c1");
    TH1D* distr = (TH1D*) c->GetPrimitive("hist")->Clone();
    distr->Scale(1/distr->Integral());
    distr->Scale(100.);
//    distr->Smooth(1);
    return distr;
}

TGraph* CamWeightNSBgraph()
{
    TGraph* gra;
    const int n = 45;
    double ang[n] = {0};
    double val[n] = {0};
    TFile* f = new TFile("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/nsb_distr_3.root");
    TCanvas* c = (TCanvas*) f->Get("c1");
    TH1D* distr = (TH1D*) c->GetPrimitive("hist")->Clone();
    distr->Scale(1/distr->Integral());
    distr->Scale(100.);
    distr->Smooth(1);
    for(int i=0; i < n; i++)
    {
        ang[i] = (double)i*2.;
        if(i > 3)
            val[i] = distr->GetBinContent(distr->FindBin(ang[i]));
    }
    ang[3] = 0.1 * ang[4]; //for better splines
    ang[2] = 0.02 * ang[4]; //for better splines
    gra = new TGraph(n, ang, val);
    return gra;
}

TH1D* CamWeightSim()
{
    TFile* f = new TFile("/Users/az/CTA/work/MThesis/mst_opt/photon_distr.root");
    TCanvas* c = (TCanvas*) f->Get("c1_n2");
    TH1D* distr = (TH1D*) c->GetPrimitive("hist")->Clone();
    distr->Scale(100.);
    return distr;
}


TGraph* SiRef()
{
    TGraph* gra = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/newrefindices/si_green.txt");
    rescaleX(gra, um / nm);
    return gra;
}

TGraph* SiO2Ref()
{
    TGraph* gra = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/newrefindices/sio2_malitson.txt");
    rescaleX(gra, um / nm);
    return gra;
}

TGraph* Si3N4Ref()
{
    TGraph* gra = new TGraph("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/newrefindices/si3n4_luke_philipp.txt");
    rescaleX(gra, um / nm);
    return gra;
}

TH1D* ledspectrum()
{
    std::string fname_str = "/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/led/spectrum" + std::to_string(intlambda) + ".txt";
    const char* fname = fname_str.c_str();
    cout << fname << endl;
    TGraph* spec = new TGraph(fname, "%lg, %lg");
    double lbin = spec->GetHistogram()->GetBinCenter(1);
    double rbin = spec->GetHistogram()->GetBinCenter(spec->GetHistogram()->GetNbinsX());
    int nbins = (int)(rbin-lbin);
    TH1D* spect_h = new TH1D("spectrum", "spectrum", nbins, lbin, rbin);
    for(int i = 0; i < nbins; i++)
    {
        double bincenter = spect_h->GetBinCenter(i);
        double value = spec->Eval(spect_h->GetBinCenter(i));
        if(value > 0.)
        {
            spect_h->Fill(bincenter, value);
            spect_h->SetBinError(i, TMath::Sqrt(value));
        }
        else
            spect_h->Fill(bincenter, 0);
    }
    spect_h->Scale(1. / spect_h->Integral());
    return spect_h;
}


//mppc geometry function
void PlaceMppc(AOpticsManager* manager, Double_t Z, const int usecoating)
{
    manager->SetVisLevel(7);
    int cellnum = (int)roundf(kActiveWidth / kCellWidth);
    //Double_t Z = pos;


    //resin material
    TGeoMaterial* mat = new TGeoMaterial("mat", 0, 0, 0);
    mat->SetTransparency(80);
    TGeoMedium* med = new TGeoMedium("med", 1, mat);

    //Si3N4 material(or not?)
    TGeoMaterial* mat1 = new TGeoMaterial("mat1", 0, 0, 0);
    mat1->SetTransparency(70);
    TGeoMedium* med1 = new TGeoMedium("med1", 1, mat1);

    //SiO2 material(or not?)
    TGeoMaterial* mat2 = new TGeoMaterial("mat2", 0, 0, 0);
    mat2->SetTransparency(60);
    TGeoMedium* med2 = new TGeoMedium("med2", 1, mat2);

    //silicon material
    TGeoMaterial* simat = new TGeoMaterial("simat", 0, 0, 0);
    simat->SetTransparency(50);
    TGeoMedium* simed = new TGeoMedium("simed", 1, simat);


//    //extract graphs
//    TFile* siX = TFile::Open("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/180323SiRefIndex/siX.root");
//    TFile* si_ = TFile::Open("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/180323SiRefIndex/si.root");
//    TFile* resinf = TFile::Open("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/resin_ref.root");
//    TCanvas *c1 = (TCanvas*)siX->Get("six_can");
//    TCanvas *c2 = (TCanvas*)si_->Get("si_can");
//    TCanvas *c3 = (TCanvas*)resinf->Get("can");
    TGraph* Sigraph = SiRef();
    TGraph* SiO2graph = SiO2Ref();
    TGraph* Si3N4graph = Si3N4Ref();
//    TGraph* resingraph = (TGraph*)c3->GetListOfPrimitives()->At(5)->Clone();
    rescaleX(Sigraph,nm);
    rescaleX(SiO2graph,nm);
    rescaleX(Si3N4graph,nm);
//    rescaleX(resingraph,nm);
  //  Si3N4graph->SetLineColor(kYellow);
  //  if(draw)
  //      c1->Draw();
  //  c2->Draw();
//    siX->Close();
//    si_->Close();
//    resinf->Close();
  //  resinf->Close();
//    delete c1;
//    delete c2;
//    delete c3;

    // Chip substrate
    TGeoPgon* pgon = new TGeoPgon(0, 360, 4, 2);
    pgon->DefineSection(0, -kResinThickness, kActiveWidth/2, kChipWidth/2);
    pgon->DefineSection(1, -kResinThickness - kChipThickness, kActiveWidth/2, kChipWidth/2);
    AObscuration* obs = new AObscuration("substrate", pgon);


    // epoxy resin
    TGeoBBox *box = new TGeoBBox("resinbox", kWidth/2, kWidth/2, kResinThickness/2);
    ALens* resin = new ALens("resin", box, med);
    resin->SetConstantRefractiveIndex(1.54);
//    resin->SetRefractiveIndex(resingraph);
    TGraph *abslen = AbslenGraph();
    resin->SetAbsorptionLength(abslen);

    // Si3N4(or not?)
    TGeoBBox *box1 = new TGeoBBox("sin3box", kWidth/2, kWidth/2, SiNThickness/2);
    ALens* sin3 = new ALens("sin3", box1, med1);
//    sin3->SetConstantRefractiveIndex(2.);
    sin3->SetRefractiveIndex(Si3N4graph);

    // SiO2(or not?)
    TGeoBBox *box2 = new TGeoBBox("sio2box", kWidth/2, kWidth/2, SiOThickness/2);
    ALens* sio2 = new ALens("sio2", box2, med2);
//    sio2->SetConstantRefractiveIndex(1.45);
    sio2->SetRefractiveIndex(SiO2graph);

    // silicon
    TGeoBBox *sibox = new TGeoBBox("sibox", kActiveWidth/2, kActiveWidth/2, kSilThickness/2);
    ALens* sil = new ALens("sil", sibox, simed);
//    sil->SetConstantRefractiveIndex(3.6);
    sil->SetRefractiveIndex(Sigraph);
    TGraph *siabslen = siAbslenGraph();
    if(absorption)
    {
        sil->SetAbsorptionLength(siabslen);
    }


    if(usecoating)
    {
        ResinThickness = kResinThickness;
        manager->GetTopVolume()->AddNode(resin, 1, new TGeoTranslation(0, 0, -kResinThickness/2 - Z));
    }
    else
        ResinThickness = 0;

    // complex MPPC cell
    int flag_grid = 0;
    TGeoBBox* cbox = new TGeoBBox("cellbox", kCellWidth/2*CellScaleFactor, kCellWidth/2*CellScaleFactor, 0.5*um);
    AObscuration* cellelectrode = new AObscuration("cellelectrode", cbox);
    TGeoBBox* grid1 = new TGeoBBox("grid1", kActiveWidth/2, kCellWidth/2*(1-CellScaleFactor), kSilThickness/2);
    AObscuration* separator1 = new AObscuration("sepa1",grid1);
    TGeoBBox* grid2 = new TGeoBBox("grid2", kCellWidth/2*(1-CellScaleFactor), kActiveWidth/2, kSilThickness/2);
    AObscuration* separator2 = new AObscuration("sepa2",grid2);
    for(int i = 0; i <= cellnum; i++)
    {
        for(int j = 0; j <= cellnum; j++)
        {
            if(i >= (cellnum/2. - 2) && i <= (cellnum/2. + 1) && j >= (cellnum/2. - 2) && j <= (cellnum/2. + 1))
            {
                Double_t x = i * kCellWidth - (kActiveWidth/2-kCellWidth/2);
                Double_t y = j * kCellWidth - (kActiveWidth/2-kCellWidth/2);
                sil->AddNode(cellelectrode,i*60+j,new TGeoTranslation(x, y, -kSilThickness+0.5*um));
            }
            else
            {
                Double_t x = i * kCellWidth - (kActiveWidth/2-kCellWidth/2);
                Double_t y = j * kCellWidth - (kActiveWidth/2-kCellWidth/2);
                Double_t z = 0;
//                if(i < cellnum && j < cellnum)
//                    sil->AddNode(cellelectrode,i*60+j,new TGeoTranslation(x, y, -kSilThickness/2 - kCellThickness/2.));
                if(flag_grid == 0)
                {

                    sil->AddNode(separator1,j,new TGeoTranslation(0,y-kCellWidth/2,z));

                }
            }

        }
        flag_grid = 1;
        sil->AddNode(separator2,i,new TGeoTranslation(i * kCellWidth - (kActiveWidth/2-kCellWidth/2)-kCellWidth/2,0,0));
    }



    //additional grid 2um

    flag_grid = 0;
    TGeoBBox* grid3 = new TGeoBBox("grid3", kActiveWidth/2, kCellWidth/2*(1-CellScaleFactor), WallThickness/2.);
    AObscuration* separator3 = new AObscuration("sepa3",grid3);
    TGeoBBox* grid4 = new TGeoBBox("grid4", kCellWidth/2*(1-CellScaleFactor), kActiveWidth/2, WallThickness/2.);
    AObscuration* separator4 = new AObscuration("sepa4",grid4);
    for(float i0 = -3.5; i0 <= 3.5; i0+=1.)
        for(float j0 = -3.5; j0 <= 3.5; j0+=1.)
        {
            Double_t X = i0*3.2*mm;
            Double_t Y = j0*3.2*mm;
            for(int i = 0; i <= cellnum; i++)
            {
                for(int j = 0; j <= cellnum; j++)
                {
                    if(i >= (cellnum/2. - 2) && i <= (cellnum/2. + 1) && j >= (cellnum/2. - 2) && j <= (cellnum/2. + 1))
                    {

                    }
                    else
                    {
                        Double_t x = i * kCellWidth - (kActiveWidth/2-kCellWidth/2);
                        Double_t y = j * kCellWidth - (kActiveWidth/2-kCellWidth/2);
                        Double_t z = - kResinThickness - SiOThickness - SiNThickness - Z + WallThickness/2. ;
                        if(flag_grid == 0)
                        {

                            manager->GetTopVolume()->AddNode(separator3,j,new TGeoTranslation(X,Y+y-kCellWidth/2,z));

                        }
                    }

                }
                flag_grid += 1;
                if (flag_grid == cellnum)
                    flag_grid = 0;
                Double_t z = - kResinThickness - SiOThickness - SiNThickness - Z + WallThickness/2.;
                manager->GetTopVolume()->AddNode(separator4,i,new TGeoTranslation(i * kCellWidth - (kActiveWidth/2-kCellWidth/2)-kCellWidth/2+X,Y,z));
            }
        }

    //cells
    TGeoBBox* cellsbox = new TGeoBBox("cellsbox", kActiveWidth/2., kActiveWidth/2., kCellThickness/2.);
    AFocalSurface* cells = new AFocalSurface("cells", cellsbox);
    //place the components and draw
    for(Int_t i = 0; i < 64; i++)//place chips
    {
        Double_t x0 = (i%8 - 3.5)*(kActiveWidth + kGap);
        Double_t y0 = (i/8 - 3.5)*(kActiveWidth + kGap);
        Double_t zz = - kResinThickness - kSilThickness - kCellThickness/2. - SiOThickness - SiNThickness - Z;
        //cells
        manager->GetTopVolume()->AddNode(cells, i+1, new TGeoTranslation(x0, y0, zz));

      Double_t z = - kResinThickness - kSilThickness/2. - SiOThickness - SiNThickness - Z;
      manager->GetTopVolume()->AddNode(sil, i + 1, new TGeoTranslation(x0, y0, z));
      manager->GetTopVolume()->AddNode(obs, i + 1, new TGeoCombiTrans(x0, y0, -kSilThickness/2 - SiOThickness - SiNThickness - Z, new TGeoRotation("", 45, 0, 0)));

    } //place chips

    //place resin & other components
//    manager->GetTopVolume()->AddNode(resin, 1, new TGeoTranslation(0, 0, -kResinThickness/2 - Z));
    manager->GetTopVolume()->AddNode(sio2, 1, new TGeoTranslation(0, 0, -kResinThickness-SiNThickness/2. - Z));
    manager->GetTopVolume()->AddNode(sin3, 1, new TGeoTranslation(0, 0, -kResinThickness- SiOThickness/2. - SiNThickness - Z));

}//mppc geometry function

//MST cone geometry function
AOpticsManager* PlaceMSTCone(int mppc = usemppc, int coat = 0, Double_t mppc_offset = 0)
{
    AOpticsManager* manager = new AOpticsManager("manager", "SC");

      // Make the world
      TGeoBBox* worldbox = new TGeoBBox("worldbox", 25*cm, 25*cm, 25*cm);
      AOpticalComponent* world = new AOpticalComponent("world", worldbox);
      manager->SetTopVolume(world);

      //define useful rotations
      TGeoRotation* rot0 = new TGeoRotation("rot0", 0, 0, 0);
      rot0->RegisterYourself();
      TGeoRotation* rot30 = new TGeoRotation("rot30", 30, 0, 0);
      rot30->RegisterYourself();
      TGeoRotation* rot60 = new TGeoRotation("rot60", 60, 0, 0);
      rot60->RegisterYourself();
      TGeoRotation* rot120 = new TGeoRotation("rot120", 120, 0, 0);
      rot120->RegisterYourself();


      //Make the cone (rot0 for edge, rot30 for vertex)
      const Double_t kRin = 24.5*mm;
      const Double_t kRout = 12.2*mm;
      const Double_t rad_glass=21.00*mm;
      const Double_t rad_cath=20.00*mm;
      const Double_t rad_PMT=15.*mm;
      const Double_t rad_window= 19.8*mm;
      const Double_t trunc = 5.52*mm;


      AGeoWinstonCone2D* coneV=0;

        coneV = new AGeoWinstonCone2D("coneV", kRin, kRout,kRin*1.102);
        TGeoPgon* pgon = new TGeoPgon("pgon", 0, 360, 6, 5);
        pgon->DefineSection(0, -coneV->GetDZ()*0.9999,     0, kRout*1.143);
        pgon->DefineSection(1, -coneV->GetDZ()*0.61,        0, kRin*0.764);
        pgon->DefineSection(2, -coneV->GetDZ()*0.,        0, kRin*0.978);
        pgon->DefineSection(3,  coneV->GetDZ()*0.625,         0, kRin*1.014);
        pgon->DefineSection(4,  coneV->GetDZ()*0.9999,     0, kRin*1.014);

        //Inner cone
        AGeoWinstonConePoly* hexV = new AGeoWinstonConePoly("hexV", 24.7754*mm, kRout, 6);

        //Cut inner cone from outer cone
//        int orientation = 0;
        TGeoCompositeShape* coneComp1a=0;
//        if (orientation == 1 ){//Edge
//          coneComp1a = new TGeoCompositeShape("coneComp1a", "pgon-hexV:rot30");
//        } else{ //Vertex
          coneComp1a = new TGeoCompositeShape("coneComp1a", "pgon:rot30-hexV:rot0");
//        }


        //Truncate Cone
        TGeoBBox* cutbox = new TGeoBBox("cutbox", 3.0*cm, 3*cm, trunc);
        TGeoTranslation* ctrunc =new TGeoTranslation("ctrunc",0,0,coneV->GetDZ()-trunc);
        ctrunc->RegisterYourself();
        TGeoCompositeShape* coneComp1 = new TGeoCompositeShape("coneComp1", "coneComp1a - cutbox:ctrunc");


        //Define reflecting surface of cone
        AMirror* coneMirror = new AMirror("coneMirror", coneComp1);

//        //Include reflectivity
//        coneMirror->SetReflectivity(graphR);

        //Include scattering uncertainty
//        int scat = 0;
//        if (scat !=0){
//          ABorderSurfaceCondition *condition = new ABorderSurfaceCondition(world, coneMirror);
//          double sigma=scat;
//          double d2r = TMath::Pi()/180.;
//          condition->SetGaussianRoughness(sigma*d2r);
//        }



        //Position cone in world
        TGeoTranslation*  ccone =new TGeoTranslation("ccone",0,0,-coneV->GetDZ()+2*trunc);


        if(mppc)
        {
            cout << mppc_offset << endl;
            PlaceMppc(manager, (2*coneV->GetDZ()-2*trunc) + mppc_offset, coat);
            cout << "mppc placed (MST cone fun)" << std::endl;
        }

        else
        {
            // Build a PMT having a curvature radius of r
            // We assume that the PMT does not have an input glass window
            TGeoSphere* sphere1 = new TGeoSphere("sphere1", fRpmt*0.9, fRpmt, 0, TMath::ASin(19*mm/fRpmt)*TMath::RadToDeg());
            AFocalSurface* fPMT = new AFocalSurface("pmt", sphere1);

            if(useAngularWeight)
            {
                TGraph* graQE = QEGraph2();
                fPMT->SetQuantumEfficiencyAngle(graQE);
            }

            Double_t Z = fRpmt + 2*(coneV->GetDZ() - trunc) - 2*mm ;
            TGeoTranslation* pmtranslation = new TGeoTranslation(0, 0, -Z);

            manager->GetTopVolume()->AddNodeOverlap(fPMT, 1, pmtranslation);
            cout << "PMT placed (MST cone fun)" << std::endl;
        }

        if(var_ref)
        {
//                    TFile f("/Users/az/CTA/work/lightguide/170411LightGuideROBAST/reflectance.root");
//                    TCanvas* can3 = (TCanvas*)f.Get("can3");
//                    gROOT->cd();
//                    TH2* h2 = (TH2*)can3->GetPrimitive("h2d_mod")->Clone("h2d_mod_clone");
            TFile f("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/2d_reflectance.root");
            TCanvas* can3 = (TCanvas*)f.Get("can");
//            TFile f("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/config/2d_ref_0.95_550.root");
//            TCanvas* can3 = (TCanvas*)f.Get("c1");
            gROOT->cd();
            TH2* h2 = (TH2*)can3->GetPrimitive("th2")->Clone("th2_clone");
            f.Close();
            h2->SetDirectory(gROOT);
            coneMirror->SetReflectance(h2);
        }
        else
            coneMirror->SetReflectance(1);

        if(cone)
        {
            world->AddNode(coneMirror, 1, ccone);
            cout << "MST cone placed" << std::endl;
        }
        manager->CloseGeometry();

        if(draw)
        {
           manager->GetTopVolume()->Draw("ogl");
        }

        return manager;
}//MST cone


Double_t rayshooter(AOpticsManager* manager, Double_t angle, Double_t wavlength, int mppc)
{
        Double_t wavelength;
        int ph_num = ph_number;
        Double_t varcounter = 0; //counter for variable wavelength
        ARayArray* array = new ARayArray();
        TRandom* rnd = new TRandom();
        Double_t phi, theta; //azimuthal & polar angles
        Double_t xs,ys; //starting position
        Double_t xd, yd, zd; //direction
        TGraph* CheGraph = CherenkovGraph();
        TGraph* nsbGraph = NSBGraph();
        TH1D* camw = CamWeightSim();
        TGraph* nsbcamw = CamWeightNSBgraph();
        TH1D* ledspec = ledspectrum();
        if(!useOdistr)
        {
            for(int i = 0; i < ph_num; i++)
            {
                if(mode >= 2 && !cone)
                {
                    if(mppc)
                    {
                        xs = rnd->Uniform(-kWidth/2 + 0.1*mm, kWidth/2 - 0.1*mm);
                        ys = rnd->Uniform(-kWidth/2 + 0.1*mm, kWidth/2 - 0.1*mm);
                    }
                    else
                    {
                        xs = rnd->Uniform(-fRpmt, fRpmt);
                        ys = rnd->Uniform(-fRpmt, fRpmt);
                        if((xs*xs + ys*ys) > fRpmt)
                        {
                            i--;
                            continue;
                        }
                    }

                }
                else
                {
                    ys = rnd->Uniform(-kSideLength,kSideLength);
                    if(ys < kOuterSideLength/2 && ys > -kOuterSideLength/2)
                    {
                        xs = rnd->Uniform(-kOuterHexSize/2, kOuterHexSize/2);
                    }
                    else if(ys < 0)
                    {
                        xs = rnd->Uniform((-ys - kOuterSideLength)*TMath::Sqrt(3), (ys + kOuterSideLength)*TMath::Sqrt(3));
                    }
                    else if(ys > 0)
                    {
                        xs = rnd->Uniform((ys - kOuterSideLength)*TMath::Sqrt(3), (-ys + kOuterSideLength)*TMath::Sqrt(3));
                    }
                    else
                    {
                        std::cerr << "ERROR in rayshooter" << std::endl;
                        continue;
                        i--;
                    }
                }
                if(!spectrum)
                    phi = 0;
                else
                    phi = rnd->Uniform(0, 360);
                if(angle < 0)
                    theta = rnd->Uniform(0, 90);
                else
                    theta = angle;

                xd = TMath::Cos(phi*TMath::DegToRad())*TMath::Sin(theta*TMath::DegToRad());
                yd = TMath::Sin(phi*TMath::DegToRad())*TMath::Sin(theta*TMath::DegToRad());
                zd = -TMath::Cos(theta*TMath::DegToRad()) ;

                if(spectrum)
                {
                    if(!filter)
                        wavelength = rnd->Uniform(285.,1000)*nm;
                    if(filter)
                        wavelength = rnd->Uniform(285.,750)*nm;
                }
                else if(mode == 1)
                {
                    wavelength = ledspec->GetRandom() * nm;
//                    wavelength = lambda;
                }
                else if(mode == 3)
                        wavelength = wavlength;

//                cout << wavelength << endl;


                ARay* ray0 = new ARay(i, wavelength, xs,ys,0, i, xd, yd, zd);
                manager->TraceNonSequential(ray0);
                Double_t var = 0; //variable for calculating photon contributions
                if(ray0->IsFocused())
                {
                        array->Add(ray0);
                        if(spectrum && cam_weight)
                            if(!nsb)
                                var = kMppc_qe*CheGraph->Eval(wavelength,0,"S")*decayfun_new(wavelength)*camw->GetBinContent(camw->FindBin(theta));
                            else
                                var = kMppc_qe*nsbGraph->Eval(wavelength,0,"")*decayfun_new(wavelength)*nsbcamw->Eval(theta, 0, "");
                        else
                            var = kMppc_qe*decayfun_new(wavelength);
                        if(mppc && filter)
                            var *= filter_gauss(wavelength);
                        varcounter += var;
                        //cout << "lambda: " << wavelength << "\tweight: " << CheGraph->Eval(wavelength,0,"S") << std::endl;
                }
                if(draw && i%500 == 0)
                {
                    TPolyLine3D* pol = ray0->MakePolyLine3D();
                    pol->SetLineWidth(2);
                    if(ray0->IsFocused())
                    {
                        pol->SetLineColor(12);
                    }
                    else
                        pol->SetLineColor(3);
                    pol->Draw();
                }

            }
            if(mode > 1 && !cone)
            {
                if(mppc)
                    cout << "square dist applied" << std::endl;
                else
                    cout << "circular dist applied" << std::endl;
            }

            else
                cout << "hexagonal dist applied" << std::endl;
        }

        else
        {
            Double_t rmax = fRin/TMath::Cos(TMath::Pi()/6);
            theta = angle;
            for(UInt_t i = 0; i < ph_num; i++)
            {
              Double_t phi = gRandom->Uniform(0, 360);

              Double_t dx =  TMath::Sin(theta*TMath::DegToRad())*TMath::Cos(phi*TMath::DegToRad());
              Double_t dy =  TMath::Sin(theta*TMath::DegToRad())*TMath::Sin(phi*TMath::DegToRad());
              Double_t dz = - TMath::Cos(theta*TMath::DegToRad());

              Double_t x = rmax;
              Double_t y = rmax;
              Double_t z = 0;
              while(x*x + y*y > rmax*rmax)
              {
                x = gRandom->Uniform(-1, 1)*rmax;
                y = gRandom->Uniform(-1, 1)*rmax;
              } // if
              ARay* ray = new ARay(0, wavelength, x, y, z, 0, dx, dy, dz);


              manager->TraceNonSequential(ray);

              if(ray->IsFocused())
              {
                      array->Add(ray);
              }
              if(draw && i%50 == 0)
              {
                  TPolyLine3D* pol = ray->MakePolyLine3D();
                  pol->SetLineWidth(2);
                  if(ray->IsFocused())
                  {
                      pol->SetLineColor(12);
                  }
                  else
                      pol->SetLineColor(3);
                  pol->Draw();
              }
            } // i

        }

        TObjArray* focused = array->GetFocused();
        int detnum = focused->GetEntries(); //focused signal

        delete rnd;


        if(usePosWeight & !mppc)
        {
          Double_t total = 0;
          TFile f("/Users/az/CTA/work/lightguide/170411LightGuideROBAST/can_average.root");
          TGraph* graph = (TGraph*)((TCanvas*)f.Get("can_average"))->GetPrimitive("anode");
          TGraph* QEgraph = hamaQEGraph();
          if(useAnodePosWeight)
          {
            Int_t n = graph->GetN();
            Double_t ymax = TMath::MaxElement(n, graph->GetY());

            for(Int_t j = 0; j < n; j++){
              Double_t y = graph->GetY()[j];
              Double_t x = graph->GetX()[j];
              graph->SetPoint(j, x, y/ymax);
            } // j
          } else
          {
            Int_t n = graph->GetN();
            for(Int_t j = 0; j < n/2; j++){
              Int_t j2 = n - j - 1;
              Double_t Y1 = graph->GetY()[j];
              Double_t Y2 = graph->GetY()[j2];
              Double_t X1 = graph->GetX()[j];
              Double_t X2 = graph->GetX()[j2];
              graph->SetPoint(j, X1, (Y1 + Y2)/2.);
              graph->SetPoint(j2, X2, (Y1 + Y2)/2.);
            } // j

            Double_t ymax = TMath::MaxElement(n, graph->GetY());

            for(Int_t j = 0; j < n; j++){
              Double_t y = graph->GetY()[j];
              Double_t x = graph->GetX()[j];
              graph->SetPoint(j, x, y/ymax);
            } // j
          }


          for(Int_t j = 0; j < detnum; j++)
          {
            ARay* ray = (ARay*)focused->At(j);
            Double_t x[3];
            ray->GetLastPoint(x);
            Double_t r = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);
            Double_t Y = fRpmt*TMath::ASin(r/fRpmt);
            Double_t weight = graph->Eval(Y/mm);
            if(gRandom->Uniform(0, 1) < weight)
            {
                if(spectrum && cam_weight)
                    if(!nsb)
                        total = total+QEgraph->Eval(ray->GetLambda(),0,"S") * CheGraph->Eval(ray->GetLambda(),0,"S") * camw->GetBinContent(camw->FindBin(theta));
                    else
                        total = total+QEgraph->Eval(ray->GetLambda(),0,"S") * nsbGraph->Eval(ray->GetLambda(),0,"") * nsbcamw->Eval(theta, 0, "");
                else
                    total = total+QEgraph->Eval(ray->GetLambda(),0,"S");
//                cout << "lambda: " << ray->GetLambda() << "\tweight: " << CheGraph->Eval(ray->GetLambda(),0,"S") << std::endl;
            } // if
          } // j
          delete array;
//          cout << "init: " << detnum << "\tweighted: " << total << std::endl;
          return total/(Double_t)ph_num*kPmt_ce/0.858514; //divide by minimum relative PDE
        }


//        output[0] = loss;
//        output[1] = rejection;

//        return (Double_t)detnum/(Double_t)ph_num;
        delete array;
        return varcounter/(Double_t)ph_num;
}


void lightguide_zenin_draw()
{

    Double_t coord[steps] = {0};
    Double_t value[steps] = {0};
    Double_t value1[steps] = {0};
    Double_t value2[steps] = {0};
    Double_t step = (finish - start)/(steps-1);
    if(mode == 0)
    {
//        gStyle->SetOptStat(0);
//        TH1D* distr = CamWeightSim();
//        distr->SetFillStyle(3006);
//        distr->SetFillColor(kBlue+1);
//        distr->SetTitle("; angle (deg); a.u.");
////        distr->Scale(1.5);

////        TH1D* nsbdistr = CamWeightNSB();
////        nsbdistr->SetFillStyle(3007);
////        nsbdistr->SetFillColor(kBlack);
////        nsbdistr->SetLineColor(kBlack);
////        nsbdistr->Scale(0.);
//        TGraph* nsbdistr = CamWeightNSBgraph();


//        TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
//        leg->AddEntry(distr, "Simulated Cherenkov ph. dist.", "f");
//        leg->AddEntry(nsbdistr, "Simulated NSB ph. dist.", "l");

//        TCanvas* c = new TCanvas();
//        nsbdistr->Draw("al*");
//        distr->Draw("same hist");
//        leg->Draw();
//        c->Update();


        //pde=========================================================
//        TGraph* pde50 = hamaPDEGraphSC();
//        TGraph* pde75 = hamaPDEGraph75SC();
//        TGraph* pde_old = OldHamaPDEGraphSC();
//        TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);

//        rescaleX(pde50, 1/nm);
//        rescaleX(pde75, 1/nm);
//        rescaleX(pde_old, 1/nm);
//        rescaleY(pde50, 100);
//        rescaleY(pde75, 100);
//        rescaleY(pde_old, 100);

//        pde50->SetTitle(";Wavelength (nm); PDE (%)");

//        pde50->SetLineColor(kBlue+1);
//        pde50->SetLineStyle(10);
//        pde50->SetLineWidth(3);
//        pde75->SetLineColor(kGreen+2);
//        pde75->SetLineStyle(9);
//        pde75->SetLineWidth(3);
//        pde_old->SetLineColor(kRed+2);
//        pde_old->SetLineStyle(0);
//        pde_old->SetLineWidth(3);

//        leg->AddEntry(pde50, "50 #mum Hamamatsu PDE data", "l");
//        leg->AddEntry(pde75, "Scaled 75 #mum Hamamatsu PDE data", "l");
//        leg->AddEntry(pde_old, "Coated 50 #mum Hamamatsu PDE data", "l");

//        TCanvas* c = new TCanvas();
//        c->SetGridx(1);
//        c->SetGridy(1);
//        pde50->Draw("al");
//        pde75->Draw("SAME L");
//        pde_old->Draw("same l");
//        leg->Draw();

//        c->Update();


//        TGraph* qe = hamaQEGraph();
//        rescaleX(qe, 1/nm);
//        rescaleY(qe, 100.);
//        qe->SetTitle(";Photon wavelength (nm); QE (%)");
//        qe->SetLineColor(kRed);
//        qe->SetLineWidth(3);
//        TCanvas* c = new TCanvas();
//        qe->Draw("al");
//        c->Update();

//        //filter-------------------
//        double wl0 = 280;
//        int wlsteps = 500;
//        double wl[wlsteps];
//        double filt[wlsteps];
//        double wl_step = 1.;
//        for(int i = 0; i < wlsteps; i++)
//        {
//            wl[i] = (wl0 + wl_step * i) * nm;
//            filt[i] = filter_gauss(wl[i]);
////            cout << filt[i] << endl;
//        }
////        cout << filter_gauss(15.) << std::endl;
//        TGraph* gg = new TGraph(wlsteps, wl, filt);
//        rescaleX(gg, 1/nm);
//        rescaleY(gg, 100);
//        gg->SetTitle("; photon wavelength (nm); transparency (%)");
////        TCanvas* cc = new TCanvas();
////        cc->SetGridx(1);
////        cc->SetGridy(1);
////        gg->Draw("al");
////        cc->Update();

//        nsb and cherenkov ================================
//        TCanvas* c1 = new TCanvas();
//        c1->SetGridx();
//        c1->SetGridy();
//        TGraph* cher = CherenkovGraph();
//        TGraph* nsb = NSBGraph();
//        TLegend *leg = new TLegend(0.5, 0.75, 0.9, 0.9);

////        TGraph* hama = hamaPDEGraph();
//        rescaleX(cher, 1/nm);
//        rescaleX(nsb, 1/nm);
//        rescaleY(nsb, 0.0035);
//        nsb->SetLineColor(kRed+1);
//        cher->SetLineColor(kViolet+2);
//        nsb->SetLineWidth(3);
//        cher->SetLineWidth(3);
//        cher->SetLineStyle(7);
//        gg->SetLineStyle(5);
//        gg->SetLineWidth(3);
//        nsb->SetTitle(";wavelength (nm); intensity (a.u.) / transmittance (%)");

//        leg->AddEntry(cher, "Cherenkov spectrum", "l");
//        leg->AddEntry(nsb, "NSB spectrum", "l");
//        leg->AddEntry(gg, "Transmittance of the filter", "l");

//        nsb->GetYaxis()->SetRangeUser(0., 1.2);

//        nsb->Draw("AL");
//        cher->Draw("same");
//        gg->Draw("same");
//        leg->Draw();

//        hama->Draw("same");

        //led spectra =================
//        TCanvas* c1 = new TCanvas();
//        TH1D* sp1 = ledspectrum(465);
////        cout << sp1->GetBinContent(sp1->FindBin(310.)) << endl;
//        TRandom* rnd = new TRandom();
//        TH1D* sp2 = new TH1D("sp", "sp", 800, 200., 1000.);
//        for(int i = 0; i < 10000; i++)
//            sp2->Fill(sp1->GetRandom());
//        sp2->Draw("e");
//        c1->Update();

//        return 0;


        //check pde at 405 nm
        AOpticsManager *co;
        co = PlaceMSTCone(1,0);
        double pde_val = rayshooter(co, 0, 405 * nm, 1);
        cout << "wl@405\t" << pde_val << endl;
        return 0;
    }

    else
    {
        if(mode == 1) //coll eff
        {
            std::string wvl = "/Users/az/CTA/work/MThesis/response/" + std::to_string(intlambda) + "_response.root";
            const char* fname = wvl.c_str();
            TFile f(fname);
            TCanvas* fcan = (TCanvas*)f.Get("c1");
            TList* list_f = (TList*)fcan->GetListOfPrimitives();
            TObject* zero = list_f->First();
            TObject* frst = list_f->After(zero);
            TGraphErrors* gre_sipm = (TGraphErrors*) list_f->After(frst);
            TGraphErrors* gre_pmt = (TGraphErrors*) list_f->After((TObject*)gre_sipm);
            Double_t exp_zeroval = gre_pmt->Eval(0., 0, "S");
            Double_t exp_zeroval_sipm = gre_sipm->Eval(0., 0, "S");
            fold(gre_sipm, sipm_fold);
            fold(gre_pmt, pmt_fold);

            TMultiGraph *mg = new TMultiGraph();
            mg->SetTitle(";Angle (deg);collection efficiency (a.u.)");
            Double_t total[3] = {0.};
//            TGraph* gr;
            TGraph* gr1;
            TGraph* gr2;
            TGraph* gr3;
            TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);

            for(int i = 0; i < steps; i++)
            {
                AOpticsManager *co;
                coord[i] = step*i;
                co = PlaceMSTCone(0,0);
                value2[i] = rayshooter(co, coord[i], lambda, 0);
                delete co;
                cout << coord[i] << "\t" << value2[i] << std::endl;
            }
            gr1 = new TGraph(steps,coord,value2);
            total[0] = gr1->Integral(0,40);
            leg->AddEntry(gr1, "PMT simulation", "l");


            gr1->SetMarkerColor(kRed+1);
            gr1->SetLineColor(kRed+1);
            gr1->SetLineWidth(3);

//            mg->Add(gr);


            leg->SetFillColor(0);

            if(usemppc)
            {
                for(int i = 0; i < steps; i++)
                {
                    AOpticsManager *co;
                    co = PlaceMSTCone(1, 0, 2*mm);
                    value1[i] = rayshooter(co, coord[i], lambda, 1);
                    delete co;
                    cout << coord[i] << "\t" << value1[i] << std::endl;
                    
                    if (value1[i] == 0 || value2[i] == 0)
                        value[i] = 0;
                    else
                        value[i] = value1[i]/value2[i];
                }
                gr2 = new TGraph(steps,coord,value1);
                total[1] = gr2->Integral(0,40);
                gr2->SetMarkerColor(kBlue+2);
                gr2->SetLineColor(kBlue+2);
                gr2->SetLineWidth(3);

                if (!spectrum)
                {
                    double scaler = exp_zeroval/gr1->Eval(0., 0, "S");
                    rescaleY(gr2, scaler);
                    rescaleY(gr1, scaler);
                }

                mg->Add(gr1);
                mg->Add(gr2);
                leg->AddEntry(gr2, "SiPM simulation", "l");
            }
            TCanvas *c1 = new TCanvas("c1","plot",800,800);
            c1->SetGridx();
            c1->SetGridy();
            mg->Draw("AC");
            mg->SetMinimum(0);
            if (!spectrum)
            {
                gre_sipm->Draw("SAMEp");
                gre_pmt->Draw("SAMEp");
                leg->AddEntry(gre_pmt, "PMT measurement", "e");
                leg->AddEntry(gre_sipm, "SiPM measurement", "e");

            }
            mg->SetMaximum(exp_zeroval_sipm * 1.3);
            leg->Draw();

            if(spectrum)
                cout << "integrated PMT: " << total[0] << "\nintegrated MPPC: "
                     << total[1]
                     << "\nrelative MPPC/PMT: " << total[1]/total[0] << std::endl;
        }

        if(mode == 3) //pde comp
        {
            TGraph* gr;
            TGraph* gr1;
            TGraph* gr2;
            TGraph* gr3;
            TGraph* gr4;
            TMultiGraph* mg = new TMultiGraph();
            mg->SetTitle(";wavelength, nm;PDE");
            TLegend *leg = new TLegend(0.65, 0.75, 0.9, 0.9);
            int steps_pde = 20;
            double coord_pde[steps_pde];
            Double_t step_l = (500.-300.)/(Double_t)(steps_pde-1);
            cout << step_l << std::endl;
            for(int i = 0; i < steps_pde; i++)
                coord_pde[i] = 0.;
//            for(int i = 0; i < steps; i++)//100µm silicone
//            {

//                AOpticsManager *co;
//                coord[i] = (300 + step_l*i)*nm;
//                co = PlaceMSTCone(1,1);

//                value[i] = rayshooter(co, inc_angle, coord[i], 1, 0);
//                delete co;

//                cout << "100s: " << coord[i] << "\t" << value[i] << std::endl;
//            }
//            gr = new TGraph(steps,coord,value);
            gr1 = hamaPDEGraph();
//            hamacut = 0;
            for(int i = 0; i <= steps_pde; i++)//100µm silicone
            {

                AOpticsManager *co;
                coord_pde[i] = (300 + step_l*i)*nm;
                co = PlaceMSTCone(1,0);

                value1[i] = rayshooter(co, inc_angle, coord_pde[i], 1);
                delete co;

                cout << "no coating: " << coord_pde[i] << "\t" << value1[i] << std::endl;
            }
            cout << "no prob here1" << std::endl;
            gr = new TGraph(steps_pde,coord_pde,value1);
            cout << gr->Eval(405 * nm) << endl;
            gr2 = hamaPDEGraphSC();
            gr3 = hamaPDEGraph75SC();
            gr4 = OldHamaPDEGraphSC();
            gr1->SetMarkerColor(kBlack);
            gr1->SetLineWidth(3);
            gr1->SetLineColor(kBlack);
            gr1->SetLineStyle(4);
            gr->SetMarkerStyle(23);
            gr->SetMarkerColor(kRed);
            gr->SetLineWidth(3);
            gr->SetLineColor(kRed);
            gr->SetLineStyle(1);
            gr->SetMarkerStyle(20);
            gr2->SetMarkerColor(kBlue);
            gr2->SetLineWidth(3);
            gr2->SetLineColor(kBlue);
            gr2->SetLineStyle(2);
            gr2->SetMarkerStyle(20);
            gr4->SetMarkerColor(kBlue+2);
            gr4->SetLineWidth(3);
            gr4->SetLineColor(kBlue+2);
            gr4->SetLineStyle(3);
            gr4->SetMarkerStyle(20);
            cout << "no prob here2" << std::endl;
            rescaleX(gr, 1/nm);
            rescaleX(gr2, 1/nm);
            rescaleX(gr4, 1/nm);
            rescaleX(gr3, 1/nm);
            mg->Add(gr);
            mg->Add(gr1);
            mg->Add(gr2);
//            mg->Add(gr3);
            mg->Add(gr4);
//            leg->AddEntry(gr, "100#mum silicone", "l");
            leg->AddEntry(gr, "Simulation", "l");
            leg->AddEntry(gr4, "50#mum coated SiPM data by Hamamatsu", "l");
            leg->AddEntry(gr2, "50#mum scaled Hamamatsu data", "l");
//            leg->AddEntry(gr3, "75#mum scaled Hamamatsu data", "l");
            cout << "no prob here3" << std::endl;

            TCanvas* c = new TCanvas();
            c->SetGridx();
            c->SetGridy();

            mg->SetMinimum(0);
            mg->SetMaximum(0.6);
            mg->Draw("AC");
            mg->GetXaxis()->SetLimits(200,1100);
            leg->Draw();
            c->Update();
            cout << "no prob here4" << std::endl;

            std::ofstream out("/Users/az/CTA/work/lightguide/sipm_lightguide_zenin/out.txt");
            std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
            std::cout.rdbuf(out.rdbuf());

            for(float i = 280.; i < 1000.; i += 1.)
            {
                double mult = gr2->Eval(i, 0, "S") / gr->Eval(i, 0, "S");
                std::cout << i << "\t" << mult << std::endl;
            }
        }//if3
        if(mode == 4) //pde comp
        {
            TGraph* gr;
            TGraph* gr1;
            TGraph* gr2;
            TGraph* gr3;
            TMultiGraph* mg = new TMultiGraph();
            mg->SetTitle(";wavelength, nm;PDE / ratio");
            TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);

            Double_t step_l = (550.-300.)/(Double_t)(steps-1);
            for(int i = 0; i < steps; i++)//100µm silicone
            {

                AOpticsManager *co;
                coord[i] = (300 + step_l*i)*nm;
                co = PlaceMSTCone(1,1);

                value[i] = rayshooter(co, inc_angle, coord[i], 1);
                delete co;

                cout << "100s: " << coord[i] << "\t" << value[i] << std::endl;
            }
            gr = new TGraph(steps,coord,value);
//            gr1 = hamaPDEGraph();
            absorption = 0;
            for(int i = 0; i < steps; i++)//100µm silicone
            {

                AOpticsManager *co;
                coord[i] = (300 + step_l*i)*nm;
                co = PlaceMSTCone(1,1);

                value1[i] = rayshooter(co, inc_angle, coord[i], 1);
                delete co;

                cout << "no absorption: " << coord[i] << "\t" << value1[i] << std::endl;
            }
            gr3 = new TGraph(steps,coord,value1);

            hamacut = 0;
            for(int i = 0; i < steps; i++)//100µm silicone
            {

                AOpticsManager *co;
                coord[i] = (300 + step_l*i)*nm;
                co = PlaceMSTCone(1,1);

                value2[i] = rayshooter(co, inc_angle, coord[i], 1);
                delete co;

                cout << "no absorption no cut: " << coord[i] << "\t" << value1[i] << std::endl;
            }
            gr1 = new TGraph(steps,coord,value2);
//            gr2 = hamaPDEGraphSC();
            gr1->SetMarkerColor(kBlack);
            gr1->SetLineWidth(3);
            gr1->SetLineColor(kBlack);
            gr1->SetLineStyle(4);
            gr->SetMarkerStyle(23);
            gr->SetMarkerColor(kRed);
            gr->SetLineWidth(3);
            gr->SetLineColor(kRed);
            gr->SetLineStyle(1);
            gr->SetMarkerStyle(20);
//            gr2->SetMarkerColor(kBlue);
//            gr2->SetLineWidth(3);
//            gr2->SetLineColor(kBlue);
//            gr2->SetLineStyle(2);
//            gr2->SetMarkerStyle(20);
            gr3->SetMarkerColor(kBlue+2);
            gr3->SetLineWidth(3);
            gr3->SetLineColor(kBlue+2);
            gr3->SetLineStyle(2);
            gr3->SetMarkerStyle(20);
            rescaleX(gr, 1/nm);
            rescaleX(gr1, 1/nm);
//            rescaleX(gr2, 1/nm);
            rescaleX(gr3, 1/nm);
            mg->Add(gr);
            mg->Add(gr1);
//            mg->Add(gr2);
            mg->Add(gr3);
            leg->AddEntry(gr, "absorbing layer + cutoff", "l");
            leg->AddEntry(gr3, "only cutoff", "l");
            leg->AddEntry(gr1, "raw curve", "l");
//            leg->AddEntry(gr2, "scaled Hamamatsu data", "l");

            TCanvas* c = new TCanvas();
            c->SetGridx();
            c->SetGridy();

            mg->SetMinimum(0);
            mg->SetMaximum(0.6);
            mg->Draw("AC");
            mg->GetXaxis()->SetLimits(200,1000);
            leg->Draw();
            c->Update();
        }//if4
        if(mode == 5) //SNR
        {
            TMultiGraph *mg = new TMultiGraph();
            mg->SetTitle(";Angle (deg);collection efficiency");
            Double_t total[3] = {0.};
//            TGraph* gr;
            TGraph* gr1;
            TGraph* gr2;
            TGraph* gr3;
            TLegend *leg = new TLegend(0.6, 0.8, 0.9, 0.9);
            //            if(useOdistr)
            //            {
            //                gr1 = plotgraph();
            //                leg->AddEntry(gr1, "PMT", "lp");
            //            }

            //            else
            //            {
                            for(int i = 0; i < steps; i++)
                            {
                                AOpticsManager *co;
                                coord[i] = step*i;
                                co = PlaceMSTCone(0,0);
                                value2[i] = rayshooter(co, coord[i], lambda, 0);
                                delete co;
                                cout << coord[i] << "\t" << value2[i] << std::endl;
                            }
                            gr1 = new TGraph(steps,coord,value2);
                            total[0] = gr1->Integral(0,29);
                            leg->AddEntry(gr1, "PMT", "l");
            //            }
//            for(int i = 0; i < steps; i++)
//            {

//                AOpticsManager *co;
//                co = PlacePMTCone(1,0);
//                value[i] = rayshooter(co, coord[i], lambda, 1);
//                delete co;

//                cout << coord[i] << "\t" << value[i] << std::endl;
//            }
//            gr = new TGraph(steps,coord,value);
//            leg->AddEntry(gr, "SiPM no coating", "lp");


            gr1->SetMarkerColor(kRed+1);
            gr1->SetLineColor(kRed+1);
            gr1->SetLineWidth(3);

//            mg->Add(gr);
            mg->Add(gr1);

            leg->SetFillColor(0);

            if(usemppc)
            {
//                //offset mppc
//                Double_t value3[steps] = {0};
//                for(int i = 0; i < steps; i++)
//                {
//                    AOpticsManager *co;
//                    co = PlaceMSTCone(1,1, 3*mm);
//                    value3[i] = rayshooter(co, coord[i], lambda, 1);
//                    delete co;
//                    cout << coord[i] << "\t" << value3[i] << std::endl;

//                }

//                gr3 = new TGraph(steps,coord,value3);
//                total[2] = gr3->Integral(0,29);
//                gr3->SetMarkerColor(kGreen+2);
//                gr3->SetLineColor(kGreen+2);
//                gr3->SetLineWidth(3);
//                mg->Add(gr3);
//                leg->AddEntry(gr3, "SiPM 100#mum coating, 3mm offset", "l");

                for(int i = 0; i < steps; i++)
                {
                    AOpticsManager *co;
                    co = PlaceMSTCone(1,0);
                    value1[i] = rayshooter(co, coord[i], lambda, 1);
                    delete co;
                    cout << coord[i] << "\t" << value1[i] << std::endl;

                    if (value1[i] == 0 || value2[i] == 0)
                        value[i] = 0;
                    else
                        value[i] = value1[i]/value2[i];
                }
                gr2 = new TGraph(steps,coord,value1);
                total[1] = gr2->Integral(0,29);
                gr2->SetMarkerColor(kBlue+2);
                gr2->SetLineColor(kBlue+2);
                gr2->SetLineWidth(3);
                mg->Add(gr2);
                leg->AddEntry(gr2, "SiPM 100#mum coating", "l");
            }
            TCanvas *c1 = new TCanvas("c1","plot",800,800);
            c1->SetGridx();
            c1->SetGridy();
            mg->Draw("AC");
            mg->SetMinimum(0);
            if(spectrum)
                mg->SetMaximum(0.25);
            else
                mg->SetMaximum(0.50);
            leg->Draw();

//            TCanvas *c2 = new TCanvas("c2", "plot2", 800, 600);
//            TGraph* gr = new TGraph(steps,coord,value);
//            gr->SetMarkerColor(kBlue+2);
//            gr->SetLineColor(kBlue+2);
//            gr->SetTitle(";Angle (deg);MPPC/PMT efficiency ratio");
//            gr->Draw("AC");

            if(spectrum)
                cout << "integrated PMT: " << total[0] << "\nintegrated MPPC: "
                     << total[1]
                     << "\nrelative MPPC/PMT (300-700): " << total[1]/total[0] << std::endl;
        }
    }

}

