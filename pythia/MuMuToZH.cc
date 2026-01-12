#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <filesystem>

using namespace Pythia8;

// nEvents jobID taskID
int main(int argc, char* argv[]) {
    Pythia pythia;


    // parse cli arguments
    int nEvents = 10000;
    int jobID = 0;
    int taskID = 0;
    if (argc > 3) {
        try {
            nEvents = std::stoi(argv[1]);
            if (nEvents <= 0) {
                std::cerr << "nEvents must be positive, got " << nEvents
                          << ". Using default 10000 instead.\n";
                nEvents = 10000;
            }
            jobID = std::stoi(argv[2]);
            if (jobID < 0) {
                std::cerr << "jobID can't be negative, got " << jobID
                          << ". Using default 0 instead.\n";
                jobID = 0;
            }
            taskID = std::stoi(argv[3]);
            if (taskID < 0) {
                std::cerr << "taskID can't be negative, got " << taskID
                          << ". Using default 0 instead.\n";
                taskID = 0;
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Could not parse integers from arguments '"
                      << argv[1] << "', '" << argv[2] << "', '" << argv[3]
	              << "'. Using defaults (nEvents = 10000, jobID = 0, taskID = 0) instead.\n";
            nEvents = 10000;
            jobID = 0;
            taskID = 0;
        }
    } else {
            std::cerr << "Using defaults (nEvents = 10000, jobID = 0, taskID = 0).\n";
    }
    std::cout << "Generating " << nEvents << " events.\n";

    // set random seed
    pythia.readString("Random:setSeed = on");
    int randSeed = jobID + taskID;
    while (randSeed > 900000000) {
        std::string s = std::to_string(randSeed);
        s.erase(0, 1);
        randSeed = std::stoi(s);
    }
    if (randSeed <= 0) randSeed = 1;  // Pythia requires positive seeds
    pythia.readString("Random:seed = " + std::to_string(randSeed));


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
    std::string outDir = "hepmc";
    std::filesystem::create_directories(outDir);
    std::string outPath = outDir + "/MuMuToZH_" + std::to_string(jobID) + "_" + std::to_string(taskID) + ".hepmc";
    HepMC3::WriterAscii writer(outPath);

    int writtenEvents = 0;
    while (writtenEvents < nEvents) {
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
	writtenEvents++;
    }

    writer.close();
    pythia.stat();
    return 0;
}
