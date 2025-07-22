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
    HepMC3::WriterAscii writer("MuMuToZH.hepmc");

    int nEvents = 10000;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if (!pythia.next()) continue;

        // Convert and write to HepMC
        HepMC3::GenEvent hepmc_evt;
        toHepMC.fill_next_event(pythia, hepmc_evt);

        // Filter particles: remove final-state particles with |eta| > 2.3
        std::vector<std::shared_ptr<HepMC3::GenParticle>> to_remove;
        for (auto p : hepmc_evt.particles())
        {
            if (p->status() != 1) continue; // Only final-state
            auto mom = p->momentum();
            double px = mom.px(), py = mom.py(), pz = mom.pz();
            double p_mag = std::sqrt(px*px + py*py + pz*pz);
            if (p_mag == std::abs(pz)) continue; // avoid div-by-zero
            double eta = 0.5 * std::log((p_mag + pz) / (p_mag - pz));
            if (std::abs(eta) > 2.3) to_remove.push_back(p);
        }

        for (auto p : to_remove)
            hepmc_evt.remove_particle(p);

        writer.write_event(hepmc_evt);

    }

    writer.close();
    pythia.stat();
    return 0;
}
