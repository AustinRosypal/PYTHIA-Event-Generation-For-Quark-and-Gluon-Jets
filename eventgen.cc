#include <iostream>
#include <vector>
#include <cmath>

#include "Pythia8/Pythia.h"

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>

using namespace Pythia8;
using std::cout;
using std::endl;
using std::vector;

// Helper function to wrap an azimuthal-angle difference into [-pi, pi].
double deltaPhi(double phi1, double phi2) {
  double dPhi = phi1 - phi2;
  while (dPhi >  M_PI) dPhi -= 2.0 * M_PI;
  while (dPhi < -M_PI) dPhi += 2.0 * M_PI;
  return dPhi;
}

// Helper function for deltaR in eta-phi space.
double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dEta = eta1 - eta2;
  double dPhi = deltaPhi(phi1, phi2);
  return std::sqrt(dEta*dEta + dPhi*dPhi);
}

int main() {

  // --------------------------------------------------------------------------
  // 1. PYTHIA setup: ep collisions at sqrt(s) = 100 GeV.
  // --------------------------------------------------------------------------
  Pythia pythia;

  pythia.readString("Beams:idA = 11");      // electron. Use -11 for positron.
  pythia.readString("Beams:idB = 2212");    // proton.
  //pythia.readString("Beams:eCM = 100.");    // center-of-mass energy in GeV.

  pythia.readString("Beams:eA = 10.");
  pythia.readString("Beams:eB = 275.");

  // Enable relevant hard processes.
  // Note: the analysis below explicitly selects subprocess code 284:
  // gamma + q -> q + g.
  pythia.readString("HardQCD:all = on");
  pythia.readString("PhaseSpace:pTHatMin = 7.");
  pythia.readString("WeakBosonExchange:all = on");
  pythia.readString("PhotonParton:all = on");
  pythia.readString("PDF:lepton2gamma = on");

  pythia.init();

  // --------------------------------------------------------------------------
  // 2. Output ROOT file and analysis configuration.
  // --------------------------------------------------------------------------
  TFile *fout = new TFile("whoa.root", "RECREATE");

  const int nEvents = 1e6;

  // SlowJet configuration.
  // power = -1 gives anti-kT-like clustering.
  // R = 0.5 is the jet radius.
  // pTmin = 3.0 GeV is the minimum reconstructed jet pT.
  const double R = 0.5;
  const double pTJetMin = 5.0;
  const double matchDRMax = R / 2.0;  // Require minDR < R/2 for truth matching.

  int num284Events = 0;
  int numQuarkJets = 0;
  int numGluonJets = 0;
  int numUnmatchedJets = 0;
  int numHardQuarks = 0;
  int numHardGluons = 0;

  double sumQuarkChargedMult = 0.0;
  double sumGluonChargedMult = 0.0;
  double sumQuarkGirth = 0.0;
  double sumGluonGirth = 0.0;

  // --------------------------------------------------------------------------
  // 3. Histograms.
  // --------------------------------------------------------------------------
  TH1F *h_quarkMult = new TH1F(
      "h_quarkMult",
      "Charged Particle Multiplicity in Quark Jets;N_{ch};Jets",
      20, 0, 20);

  TH1F *h_gluonMult = new TH1F(
      "h_gluonMult",
      "Charged Particle Multiplicity in Gluon Jets;N_{ch};Jets",
      20, 0, 20);

  TH1F *h_quarkGirth = new TH1F(
      "h_quarkGirth",
      "Quark Jet Width / Girth;g = #Sigma_{i} (p_{T,i}/p_{T,jet}) #DeltaR_{i};Jets",
      40, 0, 0.4);

  TH1F *h_gluonGirth = new TH1F(
      "h_gluonGirth",
      "Gluon Jet Width / Girth;g = #Sigma_{i} (p_{T,i}/p_{T,jet}) #DeltaR_{i};Jets",
      40, 0, 0.4);

  TH1F *h_jetpt = new TH1F(
      "h_jetpt",
      "Matched Jet p_{T} Distribution;p_{T}^{jet} [GeV];Jets",
      75, 5, 30);

  TH2F *h_quarkPhaseSpace = new TH2F(
      "h_quarkPhaseSpace",
      "Quark Jet Phase Space;#eta_{jet};#phi_{jet}",
      30, -3.5, 3.5, 30, -M_PI, M_PI);

  TH2F *h_gluonPhaseSpace = new TH2F(
      "h_gluonPhaseSpace",
      "Gluon Jet Phase Space;#eta_{jet};#phi_{jet}",
      30, -3.5, 3.5, 30, -M_PI, M_PI);

  // Struct to store outgoing hard-process partons used for truth labeling.
  struct PartonInfo {
    Vec4 p;
    int id;
  };

  // --------------------------------------------------------------------------
  // 4. Event loop.
  // --------------------------------------------------------------------------
  for (int event = 0; event < nEvents; event++) {

    if (!pythia.next()) continue;

    // Select gamma + q -> q + g events.
    int code = pythia.info.code();
    if (code != 284) continue;
    num284Events++;

    // ------------------------------------------------------------------------
    // 4a. Collect outgoing hard-process quark/gluon partons.
    //
    // IMPORTANT:
    // For many PYTHIA configurations, the cleanest hard-scattering record is
    // pythia.process, not pythia.event. In the full event record, hard-process
    // entries are often copied with negative status codes or embedded in the
    // shower history. Therefore, first use pythia.process and select entries
    // with |status| == 23, corresponding to outgoing hard-process particles.
    //
    // For subprocess 284, gamma + q -> q + g, this should usually return one
    // outgoing quark and one outgoing gluon.
    // ------------------------------------------------------------------------
    vector<PartonInfo> hardProducts;

    for (int ip = 0; ip < pythia.process.size(); ++ip) {
      int id = pythia.process[ip].id();
      int absId = std::abs(id);
      int statusAbs = std::abs(pythia.process[ip].status());

      if (statusAbs != 23) continue;

      if (absId >= 1 && absId <= 6) {
        hardProducts.push_back({pythia.process[ip].p(), id});
        numHardQuarks++;
      }
      else if (id == 21) {
        hardProducts.push_back({pythia.process[ip].p(), id});
        numHardGluons++;
      }
    }

    // Fallback: if pythia.process does not contain the hard products for some
    // reason, try the full event record using |status| == 23. This avoids the
    // all-zero problem that occurs when only status() == +23 is tested.
    if (hardProducts.empty()) {
      for (int ip = 0; ip < pythia.event.size(); ++ip) {
        int id = pythia.event[ip].id();
        int absId = std::abs(id);
        int statusAbs = std::abs(pythia.event[ip].status());

        if (statusAbs != 23) continue;

        if (absId >= 1 && absId <= 6) {
          hardProducts.push_back({pythia.event[ip].p(), id});
          numHardQuarks++;
        }
        else if (id == 21) {
          hardProducts.push_back({pythia.event[ip].p(), id});
          numHardGluons++;
        }
      }
    }

    // ------------------------------------------------------------------------
    // 4b. Cluster final-state particles into jets using SlowJet.
    // ------------------------------------------------------------------------
    SlowJet slowJet(-1, R, pTJetMin);
    slowJet.analyze(pythia.event);

    int nJets = slowJet.sizeJet();

    // If no hard products were found, do not assign quark/gluon labels. Count
    // the reconstructed jets as unmatched and move to the next event.
    if (hardProducts.empty()) {
      numUnmatchedJets += nJets;
      continue;
    }

    // 4c. Loop over reconstructed jets and match each to the nearest outgoing
    // hard-process parton. Then calculate charged multiplicity and jet girth.
    // ------------------------------------------------------------------------
    for (int j = 0; j < nJets; j++) {

      Vec4 jet = slowJet.p(j);
      vector<int> consts = slowJet.constituents(j);

      double etaJet = jet.eta();
      double phiJet = jet.phi();
      double pTJet = jet.pT();

      if (pTJet <= 0.0) continue;

      // Find nearest outgoing hard-process quark/gluon in deltaR.
      int originPDG = 0;
      double minDR = 1.0e9;

      for (const auto& parton : hardProducts) {
        double dR = deltaR(etaJet, phiJet, parton.p.eta(), parton.p.phi());

        if (dR < minDR) {
          minDR = dR;
          originPDG = parton.id;
        }
      }

      // Require the jet to be reasonably close to the matched hard parton.
      // This avoids labeling every reconstructed jet even when no good match exists.
      if (minDR >= matchDRMax) {
        numUnmatchedJets++;
        continue;
      }

      // ----------------------------------------------------------------------
      // 4d. Calculate charged particle multiplicity and jet girth.
      //
      // Charged multiplicity:
      //   N_ch = number of charged, final-state, visible particles in the jet.
      //
      // Girth / jet width:
      //   g = sum_i (pT_i / pT_jet) * DeltaR(i, jet axis),
      // using final-state visible constituents.
      // ----------------------------------------------------------------------
      int chargedMultiplicity = 0;
      double jetGirth = 0.0;

      for (int idx : consts) {
        const Particle& particle = pythia.event[idx];

        if (!particle.isFinal()) continue;
        if (!particle.isVisible()) continue;

        double dRconst = deltaR(particle.eta(), particle.phi(), etaJet, phiJet);

        // Momentum-weighted radial moment of the jet.
        jetGirth += (particle.pT() / pTJet) * dRconst;

        // Charged particle multiplicity.
        if (particle.isCharged()) chargedMultiplicity++;
      }

      // ----------------------------------------------------------------------
      // 4e. Fill quark or gluon histograms depending on matched hard parton.
      // ----------------------------------------------------------------------
      if (originPDG == 21) {
        numGluonJets++;
        sumGluonChargedMult += chargedMultiplicity;
        sumGluonGirth += jetGirth;

        h_gluonMult->Fill(chargedMultiplicity);
        h_gluonGirth->Fill(jetGirth);
        h_jetpt->Fill(pTJet);
        h_gluonPhaseSpace->Fill(etaJet, phiJet);
      }
      else if (std::abs(originPDG) <= 6) {
        numQuarkJets++;
        sumQuarkChargedMult += chargedMultiplicity;
        sumQuarkGirth += jetGirth;

        h_quarkMult->Fill(chargedMultiplicity);
        h_quarkGirth->Fill(jetGirth);
        h_jetpt->Fill(pTJet);
        h_quarkPhaseSpace->Fill(etaJet, phiJet);
      }
    }
  }

  // --------------------------------------------------------------------------
  // 5. Compute averages safely.
  // --------------------------------------------------------------------------
  double avgGluonMult = (numGluonJets > 0) ? sumGluonChargedMult / numGluonJets : 0.0;
  double avgQuarkMult = (numQuarkJets > 0) ? sumQuarkChargedMult / numQuarkJets : 0.0;
  double avgGluonGirth = (numGluonJets > 0) ? sumGluonGirth / numGluonJets : 0.0;
  double avgQuarkGirth = (numQuarkJets > 0) ? sumQuarkGirth / numQuarkJets : 0.0;

  // --------------------------------------------------------------------------
  // 6. Draw and save summary canvas.
  // --------------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Quark/Gluon Jet Analysis", 1400, 1000);
  c1->Divide(3, 3);

  c1->cd(1);
  h_quarkMult->Draw();

  c1->cd(2);
  h_gluonMult->Draw();

  c1->cd(3);
  h_quarkGirth->Draw();

  c1->cd(4);
  h_gluonGirth->Draw();

  c1->cd(5);
  h_jetpt->Draw();

  c1->cd(6);
  h_quarkPhaseSpace->Draw("colz");

  c1->cd(7);
  h_gluonPhaseSpace->Draw("colz");

  //c1->SaveAs("canvas.png");

  // --------------------------------------------------------------------------
  // 7. Write histograms to ROOT file.
  // --------------------------------------------------------------------------
  fout->cd();
  h_quarkMult->Write();
  h_gluonMult->Write();
  h_quarkGirth->Write();
  h_gluonGirth->Write();
  h_jetpt->Write();
  h_quarkPhaseSpace->Write();
  h_gluonPhaseSpace->Write();
  fout->Close();

  // --------------------------------------------------------------------------
  // 8. Print summary statistics.
  // --------------------------------------------------------------------------
  cout << endl;
  cout << "Number of generated events attempted: " << nEvents << endl;
  cout << "Number of selected subprocess-284 events: " << num284Events << endl;
  cout << "Fraction of selected subprocess-284 events: "
       << 100.0 * num284Events / nEvents << "%" << endl;
  cout << endl;

  cout << "Number of matched quark jets: " << numQuarkJets << endl;
  cout << "Number of matched gluon jets: " << numGluonJets << endl;
  cout << "Number of unmatched jets rejected by minDR < R/2: " << numUnmatchedJets << endl;
  cout << endl;

  cout << "Average charged multiplicity in quark jets: " << avgQuarkMult << endl;
  cout << "Average charged multiplicity in gluon jets: " << avgGluonMult << endl;
  cout << endl;

  cout << "Average quark jet girth: " << avgQuarkGirth << endl;
  cout << "Average gluon jet girth: " << avgGluonGirth << endl;
  cout << endl;

  cout << "Number of outgoing hard-process quarks found: " << numHardQuarks << endl;
  cout << "Number of outgoing hard-process gluons found: " << numHardGluons << endl;

  pythia.stat();

  return 0;
}

