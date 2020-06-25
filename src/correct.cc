#include <iostream>
#include <chrono>

#include <TSystem.h>
#include <GlobalConfig.h>
//#include "src/correct/CorrectionTask.h"
//#include "src/base/GlobalConfig.h"

//#include "src/base/QvectorChannelsConfig.h"
//#include "src/base/QvectorTracksConfig.h"

#include "HadesCuts.h"
#include "centrality.h"

int main(int argc, char **argv) {
  using namespace std;
  if(argc < 3){
    std::cout << "Error! Please use " << std::endl;
    std::cout << " ./correct filelist.txt corrections.root friend_filelist.txt" << std::endl;
    exit(EXIT_FAILURE);
  }
  const string fileList = argv[1];
  const string calibFileName = argv[2];

  const string event_header = "event_header";
  const string vtx_tracks = "mdc_vtx_tracks";
  const string sim_tracks = "sim_tracks";
  const string sim_event = "sim_header";

  const float beam_y = 0.74;  //TODO read from DataHeader
  // Configuration will be here
  Qn::Axis<double> pt_axis("pT", 16, 0.0, 1.6);
  Qn::Axis<double> rapidity_axis("rapidity", 15, -0.75f + beam_y, 0.75f + beam_y);

  gConfig->SetEventCuts(AnalysisTree::GetHadesEventCuts(event_header));
  AnalysisTree::Variable centrality("Centrality", event_header,
                                    {"selected_tof_rpc_hits"},
                                    [](const std::vector<double> &var){
                                      return Centrality::GetInstance()->GetCentrality5pc(var.at(0));
                                    });
  gConfig->AddEventVar(centrality);
  gConfig->AddCorrectionAxis( {"Centrality", 8, 0.0, 40.0} );

  Qn::QvectorTracksConfig un_reco("tracks_mdc", vtx_tracks,
                                         {pt_axis, rapidity_axis});
  un_reco.SetCorrectionSteps(true, true, true);
  un_reco.AddCut( {"geant_pid"}, [](double pid) { return abs(pid - 14.0) < 0.1; } );
  un_reco.SetHarmonics({1,2});
  gConfig->AddTrackQvector(un_reco);

  Qn::QvectorTracksConfig un_sim("sim_tracks", sim_tracks,
                                         {pt_axis, rapidity_axis});
  un_sim.SetCorrectionSteps(false, false, false);
  un_sim.AddCut( {"geant_pid"}, [](double pid) { return abs(pid - 14.0) < 0.1; } );
  un_sim.AddCut( {"is_primary"}, [](double flag) { return abs(flag - 1.0) < 0.1; } );
  un_sim.SetHarmonics({1,2});
  gConfig->AddTrackQvector(un_sim);

  Qn::QvectorChannelsConfig psi_rp("psi_rp", sim_event, "Ones", "reaction_plane",
                                   {1},AnalysisTree::kEventHeader);
  psi_rp.SetCorrectionSteps(false, false, false);
  gConfig->SetPsiQvector(psi_rp);

 // ***********************************************
  // first filelist should contain DataHeader

  Qn::CorrectionTask task( {fileList}, {"hades_analysis_tree"}, calibFileName);
  task.AddBranch({vtx_tracks, AnalysisTree::GetHadesTrackCuts(vtx_tracks)});
  task.AddBranch({sim_tracks});
  task.SetRecEventName("event_header");
  task.SetSimEventName("sim_header");

  std::cout << "go" << std::endl;
  auto start = std::chrono::system_clock::now();
  int n_events = task.Run();

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "elapsed time: " << elapsed_seconds.count() << " s\n";
  std::cout << ((double) n_events)*3.6 /elapsed_seconds.count() << "K events per hour\n";
  return 0;
}
