#include <TGraphErrors.h>
#include <TLegend.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TStyle.h"

struct Data {
  int N_bin;
  int N_ext;
  double Chi;
  double Chi2;
};

class Simulation {
  double k = 5.2;
  double phi = 1.8;
  double b = 0.2;
  double x_min = 0;
  double x_max = 1;
  Data data;
  int N_extract = 10000;
  int N_Bin = 100;
  std::vector<double> x, y, errY;
  std::vector<double> f_val{};
  // TF1* fitFunc;
  double k_f = 5.3;
  double phi_f = 1.8;
  double b_f = 0.2;

 public:
  // std::vector<double> f_val{};

  TF1* makeFunction() {
    TF1* f1 = new TF1("f1", "pow(cos([0]*x +[1]), 2) +[2]", x_min, x_max);
    f1->SetParameter(0, 5.2);  // gRandom->Gaus(5.2,  0.035));
    f1->SetParameter(1, 1.8);  // gRandom->Gaus(1.8,  0.03));
    f1->SetParameter(2, 0.2);  // gRandom->Gaus(0.2,  0.001));
    for (int i = 1; i <= N_Bin; ++i) {
      double fv =
          f1->Eval((x_max - x_min) / N_Bin * i - (x_max - x_min) / (2 * N_Bin));
      f_val.push_back(fv);
    }
    return f1;
  };

  TH1D* makeHisto(TF1* f1) {
    std::cout << "N_bin: " << N_Bin << ", N_extract: " << N_extract << '\n';
    TH1D* h1 = new TH1D("h1", "Histo 1", N_Bin, x_min, x_max);
    for (int i = 0; i < N_extract; ++i) {
      double x = f1->GetRandom();
      h1->Fill(x);
      h1->GetBinContent(i);
    }
    double Base = h1->GetBinWidth(1);
    double l = 0;
    for (int j = 1; j <= h1->GetNbinsX(); ++j) {
      l += h1->GetBinContent(j);
    }
    double area = f1->Integral(x_min, x_max);
    //  double xx = l / N_extract;
    double area_histo = Base * l;
    double scale_factor = area / area_histo;
    h1->Scale(scale_factor);
    // h1->Scale(xx);
    return h1;
  };

  TH1D* Generation_per_bin(TF1* f1) {
    TH1D* h2 = new TH1D("h2", "Generation per bin", N_Bin, x_min, x_max);
    for (int i = 1; i <= N_Bin; ++i) {
      double b = (x_max - x_min) / N_Bin;
      double v = f1->Eval((b)*i);
      double u = gRandom->Gaus(v, std::sqrt(v));  // std::sqrt(v));
      h2->SetBinContent(i, u);  // per asegnare al bin i in valore u
    }
    return h2;
  }

  void repete1(int N, const char* outpng = "Canva1.png") {
    std::vector<Data> vector{};
    std::vector<TH1D*> vector_histo1{};

    for (int j = 0; j < N; ++j) {
      TF1* f1 = makeFunction();
      TH1D* h1 = makeHisto(f1);
      data.N_bin = h1->GetNbinsX();
      vector_histo1.push_back(h1);
    }

    int N_Bin = vector_histo1[0]->GetNbinsX();

    std::vector<double> vector_mean(N_Bin, 0.0);
    std::vector<double> vector_dev_std(N_Bin, 0.0);

    for (int b = 1; b <= N_Bin; ++b) {
      double sum = 0.0;
      double sum_sq = 0.0;
      for (auto* h : vector_histo1) {
        double diff = h->GetBinContent(b) - vector_mean[b - 1];
        sum_sq += diff * diff;
        sum += h->GetBinContent(b);
      }
      vector_mean[b - 1] = sum / N;
      vector_dev_std[b - 1] = std::sqrt(sum_sq / (N - 1));
    }

    TCanvas* c = new TCanvas("c", "Canva 1", 800, 600);
    TH1D* fin = new TH1D("fin", "Media e Dev Std", N_Bin, x_min, x_max);

    for (int b = 1; b <= N_Bin; ++b) {
      fin->SetBinContent(b, vector_mean[b - 1]);
      fin->SetBinError(b, vector_dev_std[b - 1]);
    }

    TF1* f = makeFunction();
    double chi2 = 0;
    for (int i = 1; i <= N_Bin; ++i) {
      //  double rad = 0;

      double obs = fin->GetBinContent(i);
      double exp =
          f->Eval((x_max - x_min) / N_Bin * i - (x_max - x_min) / (2 * N_Bin));
      double rad = pow((obs - exp), 2) / exp;
      chi2 += rad;
    }
    std::cout << "CHI2_1: " << chi2 << '\n';
    fin->Draw("E1 HIST");  //"E1 HIST");
    f->Draw("SAME");
    c->SaveAs(outpng);

    std::ofstream file("data1.txt");
    if (file.is_open()) {
      file << "N_bin; mean; dev_std\n";
      for (int b = 0; b < N_Bin; ++b) {
        file << b + 1 << "; " << vector_mean[b] << "; " << vector_dev_std[b]
             << "\n";
      }
      file.close();
    } else {
      std::cerr << "Error. Unable to open file.\n";
    }

    // Pulizia memoria (importante!)
    for (auto* h : vector_histo1) delete h;
  }

  void repete2(int N, const char* outpng = "Canva2.png") {
    std::vector<Data> vector{};

    std::vector<TH1D*> vector_histo2{};

    for (int j = 0; j < N; ++j) {
      TF1* f1 = makeFunction();
      TH1D* h2 = Generation_per_bin(f1);
      data.N_bin = h2->GetNbinsX();
      vector_histo2.push_back(h2);
    }

    int N_Bin = vector_histo2[0]->GetNbinsX();

    std::vector<double> vector_mean(N_Bin, 0.0);
    std::vector<double> vector_dev_std(N_Bin, 0.0);

    for (int b = 1; b <= N_Bin; ++b) {
      double sum = 0.0;
      double sum_sq = 0.0;
      for (auto* h : vector_histo2) {
        sum += h->GetBinContent(b);
        double diff = h->GetBinContent(b) - vector_mean[b - 1];
        sum_sq += diff * diff;
      }
      vector_mean[b - 1] = sum / N;
      vector_dev_std[b - 1] = std::sqrt(sum_sq / (N - 1));
    }

    TCanvas* c = new TCanvas("c", "Canva 2", 800, 600);
    TH1D* fin = new TH1D("fin", "Media e Dev Std", N_Bin, x_min, x_max);

    for (int b = 1; b <= N_Bin; ++b) {
      fin->SetBinContent(b, vector_mean[b - 1]);
      fin->SetBinError(b, vector_dev_std[b - 1]);
    }

    TF1* f = makeFunction();
    double chi2 = 0;
    for (int i = 1; i <= N_Bin; ++i) {
      //  double rad = 0;

      double obs = fin->GetBinContent(i);
      double exp =
          f->Eval((x_max - x_min) / N_Bin * i - (x_max - x_min) / (2 * N_Bin));
      double rad = pow((obs - exp), 2) / exp;
      chi2 += rad;
    }
    std::cout << "CHI2_2: " << chi2 << '\n';
    fin->Draw("E1 HIST");  //"E1 HIST");
    f->Draw("SAME");
    c->SaveAs(outpng);

    std::ofstream file("data2.txt");
    if (file.is_open()) {
      file << "N_bin; mean; dev_std\n";
      for (int b = 0; b < N_Bin; ++b) {
        file << b + 1 << "; " << vector_mean[b] << "; " << vector_dev_std[b]
             << "\n";
      }
      file.close();
    } else {
      std::cerr << "Error. Unable to open file.\n";
    }

    // Pulizia memoria (importante!)
    for (auto* h : vector_histo2) delete h;
  }

  //------------ LAB 2---------//

  TF1* FitAnalysis() {
    TF1* fitFunc =
        new TF1("fitFunc", "pow(cos([0]*x + [1]), 2) + [2]", x_min, x_max);
    fitFunc->SetParameters(k_f, phi_f, b_f);
    fitFunc->SetParNames("k", "phi", "b");
    return fitFunc;
  }
  void LoadData(const char* filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Errore: impossibile aprire il file " << filename
                << std::endl;
      return;
    }

    std::string header;
    std::getline(file, header);  // Salta l'intestazione

    double bin, mean, sigma;
    int i = 0;
    while (file >> bin) {
      char sep;
      file >> sep >> mean >> sep >> sigma;
      x.push_back((double)i / N_Bin);  // Normalizza su [0,1]
      y.push_back(mean);
      errY.push_back(sigma);
      i++;
    }
    file.close();

    std::cout << "Letti " << x.size() << " punti da " << filename << std::endl;
  }

  void PerformFit(const char* fitname) {
    if (x.empty()) {
      std::cerr << "Nessun dato caricato!\n";
      return;
    }

    TGraphErrors* graph =
        new TGraphErrors(x.size(), x.data(), y.data(), nullptr, errY.data());
    graph->SetTitle("Fit della funzione f(x) = cos^{2}(kx + phi) + b; x; f(x)");
    graph->SetMarkerStyle(6);
    graph->SetMarkerColor(kBlue);

    TCanvas* c1 = new TCanvas("c1", "Fit", 800, 600);
    graph->Draw("AP");

    TF1* FF = FitAnalysis();
    graph->Fit(FF, "R");

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(graph, "Dati simulati", "lep");
    legend->AddEntry(FF, "Fit", "l");
    legend->Draw();

    c1->SaveAs(fitname);
  }

  void ComputeResiduals(const char* filename, const char* pngname) {
    std::ofstream file(filename);
    file << "x; residual\n";

    TCanvas* c2 = new TCanvas("c2", "Residui", 800, 600);
    TGraph* g_res = new TGraph();
    TF1* FF = FitAnalysis();
    double Res = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      double xval = x_min + (i + 0.5) * (x_max - x_min) / N_Bin;
      double f_fit = FF->Eval(xval);
      double res = y[i] - f_fit;
      g_res->SetPoint(i, x[i], res);
      file << x[i] << "; " << res << "\n";
      Res += std::abs(res);
    }
    double meanRes = Res / x.size();
    std::cout << "Residuo medio" << filename << meanRes << '\n';
    g_res->SetTitle("Residui del fit; x; (dato - fit)");
    g_res->SetMarkerStyle(104);
    g_res->SetMarkerColor(kRed);
    g_res->Draw("AP");

    TLegend* legend = new TLegend(0.55, 0.75, 0.88, 0.88);
    TString legText;
    legText.Form("Residuo medio = %.4g", meanRes);
    legend->AddEntry(g_res, legText, "p");
    legend->Draw();

    c2->SaveAs(pngname);
    file.close();
  }
  /*
    void ParameterStudy() {
      std::cout << "\n--- Studio della variazione dei parametri ---\n";
      TF1* FF = FitAnalysis();
      double k0 = FF->GetParameter(0);
      double phi0 = FF->GetParameter(1);
      double b0 = FF->GetParameter(2);

      for (double dk = -0.5; dk <= 0.5; dk += 0.05) {
        FF->SetParameter(0, k0 * (1 + dk));
        double chi2 = 0;
        for (size_t i = 0; i < x.size(); ++i) {
          double res = (y[i] - FF->Eval(x[i])) / errY[i];
          chi2 += res * res;
        }
        std::cout << "Δk/k = " << dk * 100 << "%  -> chi² = " << chi2
                  << std::endl;
      }
      std::cout << " " << '\n';
      for (double dk = -0.5; dk <= 0.5; dk += 0.05) {
        FF->SetParameter(0, phi0 * (1 + dk));
        double chi2 = 0;
        for (size_t i = 0; i < x.size(); ++i) {
          double res = (y[i] - FF->Eval(x[i])) / errY[i];
          chi2 += res * res;
        }
        std::cout << "Δphi/phi = " << dk * 100 << "%  -> chi² = " << chi2
                  << std::endl;
      }
      std::cout << " " << '\n';
      for (double dk = -0.5; dk <= 0.5; dk += 0.05) {
        FF->SetParameter(0, b0 * (1 + dk));
        double chi2 = 0;
        for (size_t i = 0; i < x.size(); ++i) {
          double res = (y[i] - FF->Eval(x[i])) / errY[i];
          chi2 += res * res;
        }
        std::cout << "Δb/b = " << dk * 100 << "%  -> chi² = " << chi2
                  << std::endl;
      }

      // Ripristina parametri
      FF->SetParameters(k0, phi0, b0);
    }*/
};

int main(int argc, char** argv) {
  Simulation sim;
  Simulation s;
  // sim.Draw("Canva1.png", "Canva1.root");
  // sim.Draw2("Canva2.png", "Canva2.root");
  sim.repete1(110);
  s.repete2(110);
  sim.FitAnalysis();
  s.FitAnalysis();
  sim.LoadData("data1.txt");
  s.LoadData("data2.txt");
  sim.PerformFit("fit1");
  s.PerformFit("fit2");
  sim.ComputeResiduals("residuals1.txt", "residuals1.png");
  s.ComputeResiduals("residuals2.txt", "residuals2.png");
  // sim.ParameterStudy();
  // s.ParameterStudy();
  std::cout << "fine " << '\n';
}
