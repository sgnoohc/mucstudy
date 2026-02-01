import pyLCIO as lcio
import argparse
import matplotlib.pyplot as plt
import numpy as np
import random
import os

COLLECTIONS = "ECalBarrelCollection,ECalEndcapCollection,HCalBarrelCollection,HCalEndcapCollection,InnerTrackerBarrelCollection,InnerTrackerEndcapCollection,OuterTrackerBarrelCollection,OuterTrackerEndcapCollection,VertexBarrelCollection,VertexEndcapCollection,YokeBarrelCollection,YokeEndcapCollection".split(',')

COLORS = "orange,orange,yellow,yellow,blue,blue,green,green,purple,purple,red,red".split(',')

TRACKERS = "InnerTrackerBarrelCollection,InnerTrackerEndcapCollection,OuterTrackerBarrelCollection,OuterTrackerEndcapCollection,VertexBarrelCollection,VertexEndcapCollection".split(',')

COLOR_MAP = dict(zip(COLLECTIONS,COLORS))

X, Y, Z = 0, 1, 2

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="Input LCIO file")
    parser.add_argument("-n", default=10, type=int,help="Max hits to dump for each collection")
    parser.add_argument("--levels", default=2, type=int,help="Number of decays to track")
    
    return parser.parse_args()

def get_collection(event, name):
    names = event.getCollectionNames()
    if name in names:
        return event.getCollection(name)
    return []

def calo_hit_from_tracked(hit, tracked_particles, min_fraction=0.0):
    total_E = hit.getEnergy()
    if total_E <= 0:
        return False

    contrib_E = 0.0

    for i in range(hit.getNMCContributions()):
        p = hit.getParticleCont(i)
        if p in tracked_particles:
            contrib_E += hit.getEnergyCont(i)

    return contrib_E / total_E >= min_fraction
    

def main():
    opts = parseArgs()
    print(f"Drawing SLCIO file {opts.i}, {opts.n} hits per each collection")
    
    reader = lcio.IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(opts.i)
    
    os.makedirs(f"./test_analysis/drawn/{opts.i[5:-6]}", exist_ok=True)
    
    for i, event in enumerate(reader):
        
        print(f"------- RUNNING EVENT {i} -------")
        
        cols = {}
        cols["MCParticle"] = get_collection(event, "MCParticle")
        
        for col in COLLECTIONS:
            cols[col] = get_collection(event, col)

            
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_xlabel('Z')
        ax.set_ylabel('X')
        ax.set_zlabel('Y')
        
        tracked_and_children = []
        particles = []
        for j, particle in enumerate(cols["MCParticle"]):
            particles.append(particle)
        
        
        goodParticle = False
        tracked_particle = None
        while not goodParticle:   
            tracked_particle = random.choice(particles)
            if tracked_particle.getGeneratorStatus() == 1 and abs(tracked_particle.getCharge()) > 1e-3:
                goodParticle = True
            

        def getTrackedChildren(level, particle):
            tracked_and_children.append(particle)
            if level > opts.levels:
                return
            children = particle.getDaughters()
            for p in children:
                getTrackedChildren(level+1, p)

        getTrackedChildren(1, tracked_particle)
        
        print(tracked_particle.getPDG())
        print(str(tracked_particle.id()) + "\n")
        
        for col in cols:
            
            if col == "MCParticle":
                continue
                
            hits_x=[]
            hits_y=[]
            hits_z=[] 
            
            for j, hit in enumerate(cols[col]):
                if col in TRACKERS:
                    if hit.getMCParticle() not in tracked_and_children:
                        continue
                else:
                    if not calo_hit_from_tracked(hit, tracked_and_children, 0.4):
                        continue
                
                hitx,hity,hitz = (hit.getPosition()[X], hit.getPosition()[Y], hit.getPosition()[Z])
                hits_x.append(hitx)
                hits_y.append(hity)
                hits_z.append(hitz)
            ax.scatter(hits_z, hits_x, hits_y, marker='o', color=COLOR_MAP[col], label=col)
        fig.savefig(f"./test_analysis/drawn/{opts.i[5:-6]}/Event {i}_oneparticle.png") 
        plt.show()
                        
if __name__ == "__main__":
    main()


