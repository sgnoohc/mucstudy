#!/usr/bin/env python3

import argparse
import math
import random
from array import array
from pyLCIO import UTIL, EVENT, IMPL, IO, IOIMPL


#_______________________________________________________________________________________________
def parse_args():
    parser = argparse.ArgumentParser(description="Sample generator steering script")

    parser.add_argument("--output", type=str, required=True, help="Path to save the generated samples")
    parser.add_argument("--nevent", type=int, default=1000, help="Number of particles to generate (default: 1000)")

    parser.add_argument("--pdg", type=int, default=13, help="PDG ID of the particle (default: 13 for muon)")
    parser.add_argument("--nparticle", type=int, default=1, help="Number of particles to generate (default: 1)")

    parser.add_argument("--pt_min", type=float, default=0.5, help="Minimum transverse momentum (GeV)")
    parser.add_argument("--pt_max", type=float, default=2.0, help="Maximum transverse momentum (GeV)")

    parser.add_argument("--eta_min", type=float, default=-2.5, help="Minimum eta")
    parser.add_argument("--eta_max", type=float, default=2.5, help="Maximum eta")

    parser.add_argument("--absphi_min", type=float, default=-math.pi, help="Minimum phi (rad)")
    parser.add_argument("--absphi_max", type=float, default=math.pi, help="Maximum phi (rad)")

    parser.add_argument("--beamspot_sigma", type=float, default=1.5, help="Beamspot Gaussian smearing sigma (in mm). Default: 1.5 mmm")


    return parser.parse_args()


#_______________________________________________________________________________________________
def get_mass(pdg_id):
    particle_masses = {
        11:   0.0005109988851472735,  # electron
        13:   0.105658,               # muon
        15:   1.77686,                # tau
        22:   0.0,                    # photon
        211:  0.13957010209560394,   # pion (π±)
        321:  0.497611,              # kaon (K±)
        2112: 0.93957,               # neutron
    }
    mass = particle_masses.get(abs(pdg_id))
    if mass is None:
        print(f"Error: Unknown PDG ID {pdg_id}", file=sys.stderr)
        sys.exit(1)
    return mass


#_______________________________________________________________________________________________
def get_charge(pdg_id):
    particle_charges = {
        11:   -1,    # electron
        13:   -1,    # muon
        15:   -1,    # tau
        22:    0,    # photon
        211:   1,    # pion+
        321:   1,    # kaon+
        2112:  0,    # neutron
    }
    abs_id = abs(pdg_id)
    charge = particle_charges.get(abs_id)
    if charge is None:
        print(f"Error: Unknown PDG ID {pdg_id}", file=sys.stderr)
        sys.exit(1)
    # Flip sign for antiparticles
    if pdg_id < 0:
        charge *= -1
    return charge


#_______________________________________________________________________________________________
def get_decay_length(pdg_id):
    decay_lengths = {
        11:   1e32,       # electron (effectively stable)
        13:   1e32,       # muon (treating as stable here)
        15:   87.03e-6,   # tau
        22:   1e32,       # photon
        211:  1e32,       # pion
        321:  15.34,      # kaon
        2112: 1e22,       # neutron
    }
    length = decay_lengths.get(abs(pdg_id))
    if length is None:
        print(f"Error: Unknown PDG ID {pdg_id}", file=sys.stderr)
        sys.exit(1)
    return length


#_______________________________________________________________________________________________
def eta_to_theta(eta):
    return 2 * math.atan(math.exp(-eta))


#_______________________________________________________________________________________________
def main():
    args = parse_args()

    print("Running with configuration:")
    for arg in vars(args):
        print(f"  {arg}: {getattr(args, arg)}")

    #--------------------------------------------

    wrt = IOIMPL.LCFactory.getInstance().createLCWriter( )
    wrt.open(args.output , EVENT.LCIO.WRITE_NEW) 
    random.seed()

    #========== particle properties ===================
    
    # particles per event
    npart = args.nparticle

    genstat  = 1

    pt_min = args.pt_min
    pt_max = args.pt_max

    theta_min = eta_to_theta(args.eta_min)
    theta_max = eta_to_theta(args.eta_max)

    pdg = args.pdg

    mass = get_mass(args.pdg)
    charge = get_charge(args.pdg)

    decayLen = get_decay_length(args.pdg)

    beamspot_sigma = args.beamspot_sigma

    #=================================================

    for j in range(0, args.nevent):

        col = IMPL.LCCollectionVec(EVENT.LCIO.MCPARTICLE) 
        evt = IMPL.LCEventImpl() 

        evt.setEventNumber(j) 

        evt.addCollection(col, "MCParticle")

        print (j, "-----------------------------")
        
        for ipart in range( 0, npart ):
        
            pt = random.uniform(pt_min, pt_max)
            # theta = random.uniform(theta_min, theta_max) 
            # phi =  random.random() * math.pi * 2.
            theta = random.uniform(140, 160) 
            phi = random.uniform(-0.2, 0.5) 

            p = pt / math.sin(theta)
            energy = math.sqrt(mass*mass + p*p) 

            if energy > 5000:
                continue
            
            px = pt * math.cos(phi)
            py = pt * math.sin(phi)
            pz = p * math.cos(theta)

            momentum  = array('f', [px, py, pz])  

            
            # --- endpoint
            epx = decayLen * math.cos(phi) * math.sin(theta) 
            epy = decayLen * math.sin(phi) * math.sin(theta)
            epz = decayLen * math.cos(theta) 

            endpoint = array('d',[epx, epy, epz])  


            # --- production vertex

            vpx = 0.
            vpy = 0.
            vpz = random.gauss(0., beamspot_sigma)

            vertex = array('d',[vpx, vpy, vpz])

            time = 0.


            # --- particle charge
            
            if ipart % 2 == 1:
                pdg = -pdg
                charge = -charge
            

            
            #--------------- create MCParticle -------------------
            
            mcp = IMPL.MCParticleImpl() 

            mcp.setGeneratorStatus(genstat) 
            mcp.setMass(mass)
            mcp.setPDG(pdg) 
            mcp.setMomentum(momentum)
            mcp.setCharge(charge) 
            mcp.setVertex(vertex)
            mcp.setTime(time)

            if decayLen < 1.e9:   # arbitrary ...
                mcp.setEndpoint(endpoint) 

            print ("  ", ipart, pdg, charge, pt, phi, theta)

            #-------------------------------------------------------

            col.addElement( mcp )

            
        wrt.writeEvent( evt ) 

    print("Generated ", pt_min, pt_max, " slice")

    wrt.close() 

if __name__ == "__main__":

    main()


