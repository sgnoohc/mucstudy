import pyLCIO as lcio
import argparse

COLLECTIONS = "ECalBarrelCollection,ECalEndcapCollection,HCalBarrelCollection,HCalEndcapCollection,InnerTrackerBarrelCollection,InnerTrackerEndcapCollection,OuterTrackerBarrelCollection,OuterTrackerEndcapCollection,VertexBarrelCollection,VertexEndcapCollection,YokeBarrelCollection,YokeEndcapCollection".split(',')

X, Y, Z = 0, 1, 2

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="Input LCIO file")
    parser.add_argument("-n", default=10, type=int,help="Max hits to dump for each collection")
    
    return parser.parse_args()

def get_collection(event, name):
    names = event.getCollectionNames()
    if name in names:
        return event.getCollection(name)
    return []


def main():
    opts = parseArgs()
    print(f"Reading SLCIO file {opts.i}, {opts.n} hits per each collection")
    
    reader = lcio.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(opts.i)
    
    #print(f"{len(reader)} event(s) found.")
    
    for i, event in enumerate(reader):
        
        print(f"------- RUNNING EVENT {i} -------")
        
        cols = {}
        cols["MCParticle"] = get_collection(event, "MCParticle")
        
        for col in COLLECTIONS:
            cols[col] = get_collection(event, col)

        for col in cols:
            print(f"{len(cols[col]):5} hits in {col}")
        print()
        
        print(f"---- Displaying first {opts.n} hits ----")
        
        for col in cols:
            if col == "MCParticle":
                continue
    
            print(f"-- {col} --\n")
            for j, hit in enumerate(cols[col]):
                if j >= opts.n:
                    break
                
                #print([method_name for method_name in dir(hit) if callable(getattr(hit, method_name))])
                
                print(f"***** HIT NO. {j+1} *****")
                if "Cal" in col or "Yoke" in col: #Calorimeters and Muon Detectors
                    print(f"Energy: {hit.getEnergy()}\nPosition: {hit.getPosition()[X]:.3e} {hit.getPosition()[Y]:.3e} {hit.getPosition()[Z]:.3e}\n")
                else: #Trackers and Vertex Detectors
                    print(f"""EDep: {hit.getEDep():.3e}
Momentum: {hit.getMomentum()[X]:.3e} {hit.getMomentum()[Y]:.3e} {hit.getMomentum()[Z]:.3e} {hit.getMomentum()[3]:.3e}\nPathLength: {hit.getPathLength():.3e}
Position: {hit.getPosition()[X]:.3e} {hit.getPosition()[Y]:.3e} {hit.getPosition()[Z]:.3e}\nQuality: {hit.getQuality()}\nTime: {hit.getTime()}\ndE/dx: {hit.getdEdx():.3e}\n""")
    
if __name__ == "__main__":
    main()


