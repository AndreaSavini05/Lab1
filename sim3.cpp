// sim3.cpp

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
  int N_extract;
  int N_Bin;

 public:
  TF1* makeFunction() {
    TF1* f1 = new TF1("f1", "pow(cos([0]*x +[1]), 2) +[2]", x_min, x_max);
    f1->SetParameter(0, k);
    f1->SetParameter(1, phi);
    f1->SetParameter(2, b);
    return f1;
  };

  TH1D* makeHisto(TF1* f1) {
    std::random_device rd;   // a seed source for the random number engine
    std::mt19937 gen(rd());  // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> ext(1000, 10000);
    std::uniform_int_distribution<> bin(70, 300);

    N_Bin = bin(gen);
    N_extract = ext(gen);
    std::cout << "N_bin: " << N_Bin << ", N_extract: " << N_extract << '\n';
    TH1D* h1 = new TH1D("h1", "Histo 1", N_Bin, x_min, x_max);
    for (int i = 0; i < N_extract; ++i) {
      double x = f1->GetRandom();
      h1->Fill(x);
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
    double ndf = h1->GetNbinsX() - f1->GetNpar();  //?
    double ki = 0;

    for (int i = 1; i <= Bin; ++i) {
      double observed = h1->GetBinContent(i);
      double Bin_mean = h1->GetBinCenter(i);
      double expected = f1->Eval(Bin_mean);

      double ki_i = pow(expected - observed, 2) / expected;
      ki += ki_i;
    }
    double chi_red = ki / ndf;

    return chi_red;
  }

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
  }

  void repete(int N) {
    std::vector<Data> vector{};
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
  sim.Draw("Canva1.png", "Canva1.root");
  sim.Draw2("Canva2.png", "Canva2.root");
  sim.repete(100);
  std::cout << "fine " << '\n';
}
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -

                                                                                    class
                                                                                    Simulation {
 
/*
Disegnare la funzione
2. Generazione ed estrazione semplice tramite estrazione random modulata secondo
la distribuzione
3. Studiare:
1. Accordo tra funzione algebrica e distribuzione generata al variare di N
eventi generati e B bins
2. Ripetere la generazione molteplici volte per stimare la l’incertezza
statistica della generazione in ogni bin e studiarne l’andamento - Incertezza da
rigenerazione.
3. Stimare l’incertezza ottenuta generando facendo fluttuare casualmente i
valori di riferimento in ogni bin con una curva Gaussiana. Incertezza da
Bin-smeering. Suggerimento procedura guidata: Generare un istogramma dalla
funzione algebrica: prendere il valore in ogni bin; farlo fluttuare in modo
Gaussiano.
4. Confrontare le incertezze ottenute.
4. Propagare le incertezze sui parametri come estrazione random attorno ad una
distribuzione Gaussiana. Considerare le seguenti incertezze per i parametri:
Delta_k +/- 2% , Delta_PHi +/-  5%, Delta_b +/- 1% . Usare i metodi 3.2 e 3.3.
*/

/*
 std::random_device seed_gen;
 std::mt19937_64 engine(seed_gen());
 std::uniform_real_distribution<> bin(10, 100);
 std::uniform_real_distribution<> ext(100, 10000);
 */
/*
TRandom3 rng(12);
double N_B = rng.Uniform(10, 100);
double N_ext = rng.Uniform(100, 1000);
int N_Bin = static_cast<double>(N_B);
int N_extract = static_cast<double>(N_ext);
*/
/*
std::random_device rd;   // a seed source for the random number engine
std::mt19937 gen(rd());  // mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<> ext(100, 1000);
std::uniform_int_distribution<> bin(10, 100);
int N_Bin = bin(gen);
int N_extract = ext(gen);*/
