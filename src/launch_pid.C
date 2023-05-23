
#include <memory>

enum class PDG {
  Proton = 2212,
  K_plus = 321,
  K_minus = -321,
  Pi_plus = 211,
  Pi_minus = -211,
  Other = -1
};

class PID {
public:
  PID(TString hist_file_name = "pid_histograms/histograms.root") {
    _hist_in = TFile::Open(hist_file_name);
    _h_proton = _hist_in->Get<TH2F>("qp_mass2_protons");
    _h_k_plus = _hist_in->Get<TH2F>("qp_mass2_K+");
    _h_k_minus = _hist_in->Get<TH2F>("qp_mass2_K-");
    _h_pi_plus = _hist_in->Get<TH2F>("qp_mass2_pi+");
    _h_pi_minus = _hist_in->Get<TH2F>("qp_mass2_pi-");
    _h_other = _hist_in->Get<TH2F>("qp_mass2_others");
    assert(_h_proton != NULL);
    assert(_h_k_plus != NULL);
    assert(_h_k_minus != NULL);
    assert(_h_pi_plus != NULL);
    assert(_h_pi_minus != NULL);
    assert(_h_pi_other != NULL);
  }
  ~PID() {
    _hist_in->Close();
    delete _h_proton;
    delete _h_k_plus;
    delete _h_k_minus;
    delete _h_pi_plus;
    delete _h_pi_minus;
  }

  PDG identify_particle(float qp, float mass2) {
    std::map<PDG, float> p_map;
    p_map.insert(
        std::make_pair(PDG::Proton, _interpolate(_h_proton, qp, mass2)));
    p_map.insert(
        std::make_pair(PDG::K_plus, _interpolate(_h_k_plus, qp, mass2)));
    p_map.insert(
        std::make_pair(PDG::K_minus, _interpolate(_h_k_minus, qp, mass2)));
    p_map.insert(
        std::make_pair(PDG::Pi_plus, _interpolate(_h_pi_plus, qp, mass2)));
    p_map.insert(
        std::make_pair(PDG::Pi_minus, _interpolate(_h_pi_minus, qp, mass2)));
    return _find_most_likely_particle(p_map);
  }

private:
  float _interpolate(TH2F *hist, float qp, float mass2) {
    auto xaxis = hist->GetXaxis();
    auto yaxis = hist->GetYaxis();
    if (qp > xaxis->GetXmin() && qp < xaxis->GetXmax() &&
        mass2 > yaxis->GetXmin() && mass2 < yaxis->GetXmax()) {
      return hist->Interpolate(qp, mass2);
    }
    return 0;
  }
  PDG _find_most_likely_particle(std::map<PDG, float> dict) {
    pair<PDG, float> entryWithMaxValue = make_pair(PDG::Other, 0);
    std::map<PDG, float>::iterator currentEntry;
    for (currentEntry = dict.begin(); currentEntry != dict.end();
         ++currentEntry) {
      if (currentEntry->second > entryWithMaxValue.second) {
        entryWithMaxValue =
            make_pair(currentEntry->first, currentEntry->second);
      }
    }
    if (entryWithMaxValue.second == 0)
      return PDG::Other;
    return entryWithMaxValue.first;
  }

  TFile *_hist_in;
  TH2F *_h_proton;
  TH2F *_h_k_plus;
  TH2F *_h_k_minus;
  TH2F *_h_pi_plus;
  TH2F *_h_pi_minus;
  TH2F *_h_other;
};

void print_mismatch(int tof_pid, int custom_pid) {
  // cout << "mismatch: tof pid: " << tof_pid << " | custom_pid: " << custom_pid
  //  << '\n';
}

void launch_pid() {

  AnalysisTree::Chain *treeIn = new AnalysisTree::Chain(
      std::vector<std::string>({"fileslist_verify.txt"}),
      std::vector<std::string>({"rTree"}));
  TFile *fileOut = TFile::Open("out/identified.root", "recreate");
  const int NEvents = treeIn->GetEntries();

  // declare branches and hook them up to the session
  auto *sim_tracks = new AnalysisTree::Particles();
  treeIn->SetBranchAddress("SimParticles.", &sim_tracks);
  auto *tof_hits = new AnalysisTree::HitDetector();
  treeIn->SetBranchAddress("TofHits.", &tof_hits);
  auto *tof_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("TofHits2SimParticles.", &tof_sim_matching);
  auto *vtx_tof_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2TofHits.", &vtx_tof_matching);
  auto *vtx_sim_matching = new AnalysisTree::Matching();
  treeIn->SetBranchAddress("VtxTracks2SimParticles.", &vtx_sim_matching);

  // declare fields to be accessed and get their id
  AnalysisTree::Configuration *treeConfig = treeIn->GetConfiguration();
  // TOF
  const int mass2 = treeConfig->GetBranchConfig("TofHits").GetFieldId("mass2");
  const int qp_tof =
      treeConfig->GetBranchConfig("TofHits").GetFieldId("qp_tof");

  // declare histograms
  TH2F hc_qp_mass2("hc_qp_mass2",
                   "correlation qp_tof mass2; sign(q)*p (GeV/c);mass^2 (GeV)^2",
                   700, -12, 12, 700, -2, 3);
  TH2F hc_qp_mass2_protons("hc_qp_mass2 protons pid",
                           "correlation qp_tof mass2 protons  pid; "
                           "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                           700, 0, 12, 700, -1.8, 3.2);
  TH2F hc_qp_mass2_pion_plus(
      "hc_qp_mass2 pi+ pid",
      "correlation qp_tof mass2 pi +  pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, 0, 10, 700, -1, 1);
  TH2F hc_qp_mass2_pion_minus(
      "hc_qp_mass2 pi-  pid",
      "correlation qp_tof mass2 pi-  pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, -10, 0, 700, -1, 1);
  TH2F hc_qp_mass2_kaon_plus(
      "hc_qp_mass2 K+  pid",
      "correlation qp_tof mass2 K+  pid; sign(q)*p (GeV/c);mass^2 (GeV)^2", 700,
      0, 10, 700, -0.8, 1);
  TH2F hc_qp_mass2_kaon_minus(
      "hc_qp_mass2 K-  pid",
      "correlation qp_tof mass2 K-  pid; sign(q)*p (GeV/c);mass^2 (GeV)^2", 700,
      -10, 0, 700, -0.8, 1);

  TH2F hc_qp_mass2_others("hc_qp_mass2 others  pid",
                          "correlation qp_tof mass2 other particles  pid; "
                          "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                          700, -12, 12, 700, -2, 5);
  TH2F hc_qp_mass2_protons_pid("hc_qp_mass2 protons sim pid",
                               "correlation qp_tof mass2 protons sim pid; "
                               "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                               700, 0, 12, 700, -1.8, 3.2);
  TH2F hc_qp_mass2_pion_plus_pid(
      "hc_qp_mass2 pi+ sim pid",
      "correlation qp_tof mass2 pi + sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, 0, 10, 700, -1, 1);
  TH2F hc_qp_mass2_pion_minus_pid(
      "hc_qp_mass2 pi- sim pid",
      "correlation qp_tof mass2 pi- sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, -10, 0, 700, -1, 1);
  TH2F hc_qp_mass2_kaon_plus_pid(
      "hc_qp_mass2 K+ sim pid",
      "correlation qp_tof mass2 K+ sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, 0, 10, 700, -0.8, 1);
  TH2F hc_qp_mass2_kaon_minus_pid(
      "hc_qp_mass2 K- sim pid",
      "correlation qp_tof mass2 K- sim pid; sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, -10, 0, 700, -0.8, 1);
  TH2F hc_qp_mass2_others_pid(
      "hc_qp_mass2 others sim pid",
      "correlation qp_tof mass2 other particles sim pid; "
      "sign(q)*p (GeV/c);mass^2 (GeV)^2",
      700, -12, 12, 700, -2, 5);

  int correctly_identified_proton = 0;
  int all_identified_proton = 0;
  int correctly_identified_k_plus = 0;
  int all_identified_k_plus = 0;
  int correctly_identified_k_minus = 0;
  int all_identified_k_minus = 0;
  int correctly_identified_pi_plus = 0;
  int all_identified_pi_plus = 0;
  int correctly_identified_pi_minus = 0;
  int all_identified_pi_minus = 0;
  int correctly_identified_other = 0;
  int all_identified_other = 0;

  // fill histograms
  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

    for (const auto &tof_hit : *tof_hits) {
      const float tof_mass2 = tof_hit.GetField<float>(mass2);
      const float tof_qp_tof = tof_hit.GetField<float>(qp_tof);

      hc_qp_mass2.Fill(tof_qp_tof, tof_mass2);

      // get real pid
      const int tof2sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      const int vtx_id = vtx_tof_matching->GetMatch(tof_hit.GetId());
      const int vtx2sim_id = vtx_sim_matching->GetMatch(vtx_id);

      if (tof2sim_id < 0 || vtx2sim_id < 0)
        continue;

      // assert that pdg from vtx is the same as the one from tof,
      // we get less data this way but it makes more sense this way and the
      // identification is more exact.
      const int tof_pdg = sim_tracks->GetChannel(tof2sim_id).GetPid();
      const int vtx_pdg = sim_tracks->GetChannel(vtx2sim_id).GetPid();
      if (tof_pdg != vtx_pdg)
        continue;

      // get my pid
      static PID pid;
      int pid_pdg = (int)pid.identify_particle(tof_qp_tof, tof_mass2);

      switch (pid_pdg) {
      case 2212: // protons
        hc_qp_mass2_protons.Fill(tof_qp_tof, tof_mass2);
        all_identified_proton++;
        break;
      case 321:
        hc_qp_mass2_kaon_plus.Fill(tof_qp_tof, tof_mass2);
        all_identified_k_plus++;
        break;
      case -321:
        hc_qp_mass2_kaon_minus.Fill(tof_qp_tof, tof_mass2);
        all_identified_k_minus++;
        break;
      case 211:
        hc_qp_mass2_pion_plus.Fill(tof_qp_tof, tof_mass2);
        all_identified_pi_plus++;
        break;
      case -211:
        hc_qp_mass2_pion_minus.Fill(tof_qp_tof, tof_mass2);
        all_identified_pi_minus++;
        break;
      default:
        hc_qp_mass2_others.Fill(tof_qp_tof, tof_mass2);
        all_identified_other++;
      }

      switch (tof_pdg) {
      case 2212: // protons
        if (tof_qp_tof < 2 && tof_mass2 < 0.6)
          goto default;
        if (tof_qp_tof < 4 && tof_mass2 < 0.4)
          continue;
        if (tof_qp_tof < 6 && tof_mass2 < 0.2)
          continue;
        if (tof_qp_tof < 4 && tof_mass2 > 1.4)
          continue;
        if (tof_qp_tof < 6 && tof_mass2 > 1.6)
          continue;
        if (tof_qp_tof < 8 && tof_mass2 > 2.2)
          continue;
        hc_qp_mass2_protons_pid.Fill(tof_qp_tof, tof_mass2);
        if (tof_pdg == pid_pdg)
          correctly_identified_proton++;
        else
          print_mismatch(tof_pdg, pid_pdg);
        break;
      case 321:
        if (tof_qp_tof > 2 && tof_mass2 > 0.14)
          hc_qp_mass2_kaon_plus_pid.Fill(tof_qp_tof, tof_mass2);
        if (tof_pdg == pid_pdg)
          correctly_identified_k_plus++;
        else
          print_mismatch(tof_pdg, pid_pdg);
        break;
      case -321:
        hc_qp_mass2_kaon_minus_pid.Fill(tof_qp_tof, tof_mass2);
        if (tof_pdg == pid_pdg)
          correctly_identified_k_minus++;
        else
          print_mismatch(tof_pdg, pid_pdg);
        break;
      case 211:
        if (tof_qp_tof < 5 && tof_mass2 > 0.4)
          continue;
        if (tof_qp_tof < 3 && tof_mass2 > 0.2)
          continue;
        if (tof_qp_tof < 1 && tof_mass2 > 0.1)
          continue;
        hc_qp_mass2_pion_plus_pid.Fill(tof_qp_tof, tof_mass2);
        if (tof_pdg == pid_pdg)
          correctly_identified_pi_plus++;
        else
          print_mismatch(tof_pdg, pid_pdg);
        break;
      case -211:
        if (tof_qp_tof > -5 && tof_mass2 > 0.4)
          continue;
        if (tof_qp_tof > -3 && tof_mass2 > 0.2)
          continue;
        if (tof_qp_tof > -1 && tof_mass2 > 0.1)
          continue;
        hc_qp_mass2_pion_minus_pid.Fill(tof_qp_tof, tof_mass2);
        if (tof_pdg == pid_pdg)
          correctly_identified_pi_minus++;
        else
          print_mismatch(tof_pdg, pid_pdg);
        break;
      default:
        hc_qp_mass2_others_pid.Fill(tof_qp_tof, tof_mass2);
        if (pid_pdg == -1)
          correctly_identified_other++;
        else
          print_mismatch(tof_pdg, pid_pdg);
      }
    }
  }
  cout << "all identified protons: " << all_identified_proton << "\n";
  cout << "all identified K+: " << all_identified_k_plus << "\n";
  cout << "all identified K-: " << all_identified_k_minus << "\n";
  cout << "all identified Pi+: " << all_identified_pi_plus << "\n";
  cout << "all identified Pi-: " << all_identified_pi_minus << "\n";
  cout << "all identified others: " << all_identified_other << "\n\n";

  cout << "correctly identified protons: " << correctly_identified_proton
       << "\n";
  cout << "correctly identified K+: " << correctly_identified_k_plus << "\n";
  cout << "correctly identified K-: " << correctly_identified_k_minus << "\n";
  cout << "correctly identified Pi+: " << correctly_identified_pi_plus << "\n";
  cout << "correctly identified Pi-: " << correctly_identified_pi_minus << "\n";
  cout << "correctly identified others: " << correctly_identified_other
       << "\n\n";

  cout << "proton purity: "
       << (float)correctly_identified_proton / (float)all_identified_proton
       << "\n";
  cout << "K+ purity: "
       << (float)correctly_identified_k_plus / (float)all_identified_k_plus
       << "\n";
  cout << "K- purity: "
       << (float)correctly_identified_k_minus / (float)all_identified_k_minus
       << "\n";
  cout << "Pi+ purity: "
       << (float)correctly_identified_pi_plus / (float)all_identified_pi_plus
       << "\n";
  cout << "Pi- purity: "
       << (float)correctly_identified_pi_minus / (float)all_identified_pi_minus
       << "\n";
  cout << "Other purity: "
       << (float)correctly_identified_other / (float)all_identified_other
       << "\n";

  // Normalize histograms
  //   hc_qp_mass2.Scale(1. / hc_qp_mass2.Integral());
  //   hc_qp_mass2_protons.Scale(1. / hc_qp_mass2_protons.Integral());
  //   hc_qp_mass2_kaon_plus.Scale(1. / hc_qp_mass2_kaon_minus.Integral());
  //   hc_qp_mass2_kaon_minus.Scale(1. / hc_qp_mass2_kaon_minus.Integral());
  //   hc_qp_mass2_pion_plus.Scale(1. / hc_qp_mass2_pion_plus.Integral());
  //   hc_qp_mass2_pion_minus.Scale(1. / hc_qp_mass2_pion_minus.Integral());
  //   hc_qp_mass2_others.Scale(1. / hc_qp_mass2_others.Integral());

  // write to histograms
  fileOut->cd();

  hc_qp_mass2.Write();

  hc_qp_mass2_protons.Write();
  hc_qp_mass2_kaon_plus.Write();
  hc_qp_mass2_kaon_minus.Write();
  hc_qp_mass2_pion_plus.Write();
  hc_qp_mass2_pion_minus.Write();
  hc_qp_mass2_others.Write();

  hc_qp_mass2_protons_pid.Write();
  hc_qp_mass2_kaon_plus_pid.Write();
  hc_qp_mass2_kaon_minus_pid.Write();
  hc_qp_mass2_pion_plus_pid.Write();
  hc_qp_mass2_pion_minus_pid.Write();
  hc_qp_mass2_others_pid.Write();

  fileOut->Close();
}