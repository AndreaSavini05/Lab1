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

 public:
  TF1* makeFunction() {
    TF1* f1 = new TF1("f1", "pow(cos([0]*x +[1]), 2) +[2]", x_min, x_max);
    f1->SetParameter(0, k);
    f1->SetParameter(1, phi);
    f1->SetParameter(2, b);

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
    double area_histo = Base * l;
    double scale_factor = area / area_histo;
    h1->Scale(scale_factor);
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

  double ChiSquare(TF1* f1, TH1D* h1) {
    int Bin = h1->GetNbinsX();
    // double ndf = h1->GetNbinsX() - f1->GetNpar();  //?
    double ki = 0;

    for (int i = 1; i <= Bin; ++i) {
      double observed = h1->GetBinContent(i);
      double Bin_mean = h1->GetBinCenter(i);
      double expected = f1->Eval(Bin_mean);

      double ki_i = pow(expected - observed, 2);  // / expected;
      ki += ki_i;
    }
    double chi_red = ki / (Bin * N_extract);

    return chi_red;
  } /*
    void Draw(const char* outpng = "Canva1.png",
               const char* outroot = "Canva1.root") {
       TCanvas* c = new TCanvas("c", "Canva 1", 800, 600);

       TF1* f1 = makeFunction();
       f1->Draw();

       TH1D* h1 = makeHisto(f1);  //, N_Bin, N_extract);
       TH1D* h2 = Generation_per_bin(f1);
       double chi = ChiSquare(f1, h1);
       double chi2 = ChiSquare(f1, h2);
       h1->Draw("HIST SAME");
       c->Update();
       c->SaveAs(outpng);

       TFile fout(outroot,
                  "RECREATE");  // recreate serve per cancellare contenuto
                                // gia presente per scriverne di nuovo
       f1->Write("f_analityc");
       c->Write("Canvas");
       fout.Close();

       std::cout << "Reducted chi square: " << chi << ", " << chi2 << '\n';
       std::cout << "Immagine creata come canva 1; file salvato in canva1.png"
                 << '\n';
     }

     void Draw2(const char* out_png = "canva2.png",
                const char* out_root = "Canva2.root") {
       TCanvas* g = new TCanvas("g", "Canva 2", 800, 600);

       TF1* f1 = makeFunction();
       f1->Draw();
       TH1D* h2 = Generation_per_bin(f1);
       double chi2 = ChiSquare(f1, h2);
       h2->Draw("HIST SAME");
       g->Update();
       g->SaveAs(out_png);

       TFile fout(out_root,
                  "RECREATE");  // recreate serve per cancellare contenuto
                                // gia presente per scriverne di nuovo
       f1->Write("f_analityc");
       g->Write("Canvas");
       fout.Close();

       std::cout << "Reducted chi square:" << chi2 << '\n';
       std::cout << "Immagine creata come canva 1; file salvato in canva2.png"
                 << '\n';
     }*/

  void repete(int N, const char* outpng = "Canva1.png") {
    std::vector<Data> vector{};
    std::vector<double> vector_mean{};
    std::vector<double> vector_dev_std{};
    for (int j = 1; j <= N; ++j) {
      TF1* f1 = makeFunction();
      TH1D* h1 = makeHisto(f1);
      TH1D* h2 = Generation_per_bin(f1);
      double chi = ChiSquare(f1, h1);
      double chi2 = ChiSquare(f1, h2);
      data.N_bin = h1->GetNbinsX();
      data.N_ext = N_extract;
      data.Chi = chi;
      data.Chi2 = chi2;
      vector.push_back(data);
    }

    for (int i = 1; i <= N_Bin; ++i) {
      double l = 0;
      double t = 0;
      double cb = h1->GetBinCentr(i);
      double f_v = f1->Eval(cb);
      for (int j = 1; j <= 100; ++j) {
        double bin_val = h1->GetBinContent(i);
        l += bin_val;
        double mean = l / 100;  // numero istogrammi

        t += pow(f_v - bin_val, 2);
        double dev_std = std::sqrt(t);
        vector_mean.push_back(mean);
        vector_dev.push_back(dev_std);
      }
    }
    TCanvas* c = new TCanvas("c", "Canva 1", 800, 600);
    TH1D* fin = new TH1D("fin", "FInale", N_Bin, x_min, x_max);
    for (int i = 1; i <= N_Bin; ++i) {
      fin->SetBinContent(i, vector_mean[i - 1]);
      fin->SetBinError(i, vector_dev_std[i - 1]);
    }
    f1->Draw();

    fin->Draw("HIST SAME", "E");
    c->Update();
    c->SaveAs(outpng);

    std::ofstream file("data.txt");
    if (file.is_open()) {
      file << " N bin " << " " << "N_extract" << " " << "chi red1" << "chiRed_2"
           << '\n';
      std::for_each(vector.begin(), vector.end(), [&](const auto& data) {
        file << data.N_bin << "; " << data.N_ext << "; " << data.Chi << "; "
             << data.Chi2 << '\n';
      });
    } else {
      std::cerr << "Error. Unable to open the file" << "\n";
    }
    file.close();
  }
};

int main(int argc, char** argv) {
  Simulation sim;
  // sim.Draw("Canva1.png", "Canva1.root");
  // sim.Draw2("Canva2.png", "Canva2.root");
  sim.repete(100);
  std::cout << "fine " << '\n';
}

/*
  for (int i = 1; i <= 100; ++i){
    for (int j = 1; j <= N_Bin; ++j){

     vector_val(i).push_back(j);
    }
  }*/
