void generate_histograms() {
  AnalysisTree::Chain *treeIn =
      new AnalysisTree::Chain(std::vector<std::string>({"fileslist.txt"}),
                              std::vector<std::string>({"rTree"}));
  TFile *fileOut = TFile::Open("pid_histograms/histograms.root", "recreate");
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
  TH2F hc_qp_mass2("qp_mass2",
                   "correlation qp_tof mass2; sign(q)*p (GeV/c);mass^2 (GeV)^2",
                   700, -12, 12, 700, -2, 3);
  TH2F hc_qp_mass2_protons("qp_mass2_protons",
                           "correlation qp_tof mass2 protons; "
                           "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                           700, 0, 12, 700, -1.8, 3.2);
  TH2F hc_qp_mass2_pion_plus(
      "qp_mass2_pi+",
      "correlation qp_tof mass2 pi +; sign(q)*p (GeV/c);mass^2 (GeV)^2", 700, 0,
      10, 700, -1, 1);
  TH2F hc_qp_mass2_pion_minus(
      "qp_mass2_pi-",
      "correlation qp_tof mass2 pi-; sign(q)*p (GeV/c);mass^2 (GeV)^2", 700,
      -10, 0, 700, -1, 1);
  TH2F hc_qp_mass2_kaon_plus(
      "qp_mass2_K+",
      "correlation qp_tof mass2 K+; sign(q)*p (GeV/c);mass^2 (GeV)^2", 700, 0,
      10, 700, -0.8, 1);
  TH2F hc_qp_mass2_kaon_minus(
      "qp_mass2_K-",
      "correlation qp_tof mass2 K-; sign(q)*p (GeV/c);mass^2 (GeV)^2", 700, -10,
      0, 700, -0.8, 1);
  TH2F hc_qp_mass2_others("qp_mass2_others",
                          "correlation qp_tof mass2 other particles; "
                          "sign(q)*p (GeV/c);mass^2 (GeV)^2",
                          700, -12, 12, 700, -2, 5);

  // fill histograms
  for (int i = 0; i < NEvents; i++) {
    treeIn->GetEntry(i);

    for (const auto &tof_hit : *tof_hits) {
      const float tof_mass2 = tof_hit.GetField<float>(mass2);
      const float tof_qp_tof = tof_hit.GetField<float>(qp_tof);

      hc_qp_mass2.Fill(tof_qp_tof, tof_mass2);

      const int tof2sim_id = tof_sim_matching->GetMatch(tof_hit.GetId());
      const int vtx_id = vtx_tof_matching->GetMatch(tof_hit.GetId());
      const int vtx2sim_id = vtx_sim_matching->GetMatch(vtx_id);

      if (tof2sim_id < 0 || vtx2sim_id < 0)
        continue;
      if (tof2sim_id != vtx2sim_id)
        continue;
      // match track to proton
      const int tof_pdg = sim_tracks->GetChannel(tof2sim_id).GetPid();
      const int vtx_pdg = sim_tracks->GetChannel(vtx2sim_id).GetPid();

      if (tof_pdg != vtx_pdg)
        continue;

      switch (tof_pdg) {
      case 2212: // protons
        if (tof_qp_tof < 2 && tof_mass2 < 0.6)
          continue;
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
        hc_qp_mass2_protons.Fill(tof_qp_tof, tof_mass2);
        break;
      case 321:
        if (tof_qp_tof < 2 && tof_mass2 < 0.14)
          continue;
        hc_qp_mass2_kaon_plus.Fill(tof_qp_tof, tof_mass2);
        break;
      case -321:
        hc_qp_mass2_kaon_minus.Fill(tof_qp_tof, tof_mass2);
        break;
      case 211:
        if (tof_qp_tof < 5 && tof_mass2 > 0.4)
          continue;
        if (tof_qp_tof < 3 && tof_mass2 > 0.2)
          continue;
        if (tof_qp_tof < 1 && tof_mass2 > 0.1)
          continue;
        hc_qp_mass2_pion_plus.Fill(tof_qp_tof, tof_mass2);
        break;
      case -211:
        if (tof_qp_tof > -5 && tof_mass2 > 0.4)
          continue;
        if (tof_qp_tof > -3 && tof_mass2 > 0.2)
          continue;
        if (tof_qp_tof > -1 && tof_mass2 > 0.1)
          continue;
        hc_qp_mass2_pion_minus.Fill(tof_qp_tof, tof_mass2);
        break;
      default:
        hc_qp_mass2_others.Fill(tof_qp_tof, tof_mass2);
      }
    }
  }

  // write to histograms
  // hc_qp_mass2.Scale(1. / hc_qp_mass2.Integral());
  // hc_qp_mass2_protons.Scale(1. / hc_qp_mass2_protons.Integral());
  // hc_qp_mass2_kaon_plus.Scale(1. / hc_qp_mass2_kaon_minus.Integral());
  // hc_qp_mass2_kaon_minus.Scale(1. / hc_qp_mass2_kaon_minus.Integral());
  // hc_qp_mass2_pion_plus.Scale(1. / hc_qp_mass2_pion_plus.Integral());
  // hc_qp_mass2_pion_minus.Scale(1. / hc_qp_mass2_pion_minus.Integral());
  // hc_qp_mass2_others.Scale(1. / hc_qp_mass2_others.Integral());

  hc_qp_mass2.Write();
  hc_qp_mass2_protons.Write();
  hc_qp_mass2_kaon_plus.Write();
  hc_qp_mass2_kaon_minus.Write();
  hc_qp_mass2_pion_plus.Write();
  hc_qp_mass2_pion_minus.Write();
  hc_qp_mass2_others.Write();

  fileOut->Close();
}