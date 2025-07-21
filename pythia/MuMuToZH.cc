#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

using namespace Pythia8;

int main() {
    Pythia pythia;

    // Set μ⁺μ⁻ beams
    pythia.readString("Beams:idA = -13");
    pythia.readString("Beams:idB = 13");
    pythia.readString("Beams:eCM = 10000.");  // sqrt(s)s = 10 TeV

    // Enable HZ production
    pythia.readString("HiggsSM:ffbar2HZ = on");

    // Optional: visible decays
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5 11 13"); // hadronic & leptonic Z decays
    pythia.readString("25:onMode = off");
    pythia.readString("25:onIfAny = 1 2 3 4 5 11 13 15 21 22"); // visible Higgs decays

    // Initialize
    pythia.init();

    // HepMC3 output
    HepMC3::Pythia8ToHepMC3 toHepMC;
    HepMC3::WriterAscii writer("MuMuToZH_ZPt3to4_HPt3to4.hepmc");

    int nEvents = 10000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        // bool keep_event = false;
        // double pT_Z = -1.0;
        // double pT_H = -1.0;

        // for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {
        //     const auto& p = pythia.event[iPart];
        //     if (p.id() == 23) pT_Z = p.pT();
        //     if (p.id() == 25) pT_H = p.pT();
        // }

        // // Require both to be in 3–4 TeV window
        // if (pT_Z > 3000. && pT_Z < 4000. &&
        //     pT_H > 3000. && pT_H < 4000.) {
        //     keep_event = true;
        // }

        // if (!keep_event) continue;

        // Convert and write to HepMC
        HepMC3::GenEvent hepmc_evt;
        toHepMC.fill_next_event(pythia, hepmc_evt);
        writer.write_event(hepmc_evt);
    }

    writer.close();
    pythia.stat();
    return 0;
}
