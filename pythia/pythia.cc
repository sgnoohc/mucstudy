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
    std::string outDir = "";
    if (argc < 6) {
        std::cout << "Usage: " << argv[0] << " <pythia_config_file> <nEvents> <jobID> <taskID> <outDir>" << std::endl;
        exit(-1);
    }
    if (argc >= 6) {
        try {
            pythia.readFile(argv[1]);

            nEvents = std::stoi(argv[2]);
            if (nEvents <= 0) {
                std::cerr << "nEvents must be positive, got " << nEvents
                          << ". Using default 10000 instead.\n";
                nEvents = 10000;
            }
            jobID = std::stoi(argv[3]);
            if (jobID < 0) {
                std::cerr << "jobID can't be negative, got " << jobID
                          << ". Using default 0 instead.\n";
                jobID = 0;
            }
            taskID = std::stoi(argv[4]);
            if (taskID < 0) {
                std::cerr << "taskID can't be negative, got " << taskID
                          << ". Using default 0 instead.\n";
                taskID = 0;
            }
	    outDir = argv[5];
        }
        catch (const std::exception& e) {
            std::cout << "Usage: " << argv[0] << " <pythia_config_file> <nEvents> <jobID> <taskID> <outDir>" << std::endl;
            exit(-1);
        }
    } else {
        std::cout << "Usage: " << argv[0] << " <pythia_config_file> <nEvents> <jobID> <taskID> <outDir>" << std::endl;
        exit(-1);
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

    pythia.readFile("pythia.conf");

    // Initialize
    pythia.init();

    // HepMC3 output
    HepMC3::Pythia8ToHepMC3 toHepMC;
    std::filesystem::create_directories(outDir);
    std::filesystem::path cfgPath(argv[1]);
    std::string cfgTag = cfgPath.stem().string();  // filename without .conf
    std::string outPath = outDir + "/" + cfgTag + "_" + std::to_string(jobID) + "_" + std::to_string(taskID) + ".hepmc";
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
