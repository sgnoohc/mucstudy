import pyLCIO as lcio
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

COLLECTIONS = "ECalBarrelCollection,ECalEndcapCollection,HCalBarrelCollection,HCalEndcapCollection,InnerTrackerBarrelCollection,InnerTrackerEndcapCollection,OuterTrackerBarrelCollection,OuterTrackerEndcapCollection,VertexBarrelCollection,VertexEndcapCollection,YokeBarrelCollection,YokeEndcapCollection".split(',')

COLORS = "orange,orange,yellow,yellow,blue,blue,green,green,purple,purple,red,red".split(',')

COLOR_MAP = dict(zip(COLLECTIONS,COLORS))

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
        ax.set_title(f"Hits on Event {i}")
        ax.set_xlabel('Z')
        ax.set_ylabel('X')
        ax.set_zlabel('Y')
        
        for col in cols:
            if col == "MCParticle":
                continue

            hits_x=[]
            hits_y=[]
            hits_z=[]
            
            for j, hit in enumerate(cols[col]):
                if j >= opts.n:
                    break

                hitx,hity,hitz = (hit.getPosition()[X], hit.getPosition()[Y], hit.getPosition()[Z])
                hits_x.append(hitx)
                hits_y.append(hity)
                hits_z.append(hitz)
            ax.scatter(hits_z, hits_x, hits_y, marker='o', color=COLOR_MAP[col], label=col)
        fig.savefig(f"./test_analysis/drawn/{opts.i[5:-6]}/Event {i}.png") 
        plt.show()
                        
            
             
                
                

if __name__ == "__main__":
    main()


