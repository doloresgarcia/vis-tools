#! /usr/bin/env python3
from podio import root_io
from particle import Particle
import sys

def get_sim_status_description():
    return (
        "simulator status bits: [sbvtcls] "
        " s: created in simulation"
        " b: backscatter"
        " v: vertex is not endpoint of parent"
        " t: decayed in tracker"
        " c: decayed in calorimeter"
        " l: has left detector"
        " s: stopped"
        " o: overlay\n"
    )

def get_sim_status(mc):
    status = list("[    0   ]")

    if mc.getSimulatorStatus() == 0:
        return "".join(status)

    flags = [
        (1, mc.isCreatedInSimulation(), 's'),
        (2, mc.isBackscatter(), 'b'),
        (3, mc.vertexIsNotEndpointOfParent(), 'v'),
        (4, mc.isDecayedInTracker(), 't'),
        (5, mc.isDecayedInCalorimeter(), 'c'),
        (6, mc.hasLeftDetector(), 'l'),
        (7, mc.isStopped(), 's'),
        (8, mc.isOverlay(), 'o'),
    ]

    for idx, condition, char in flags:
        status[idx] = char if condition else ' '

    return "".join(status)

def mc_print_header():
    header = (
        f'{"index":^5}|'
        f'{"PDG":^10}|'
        f'{"name":^15}|'
        f'{"px, py, pz (GeV)":^25}|'
        f'{"px_ep, py_ep, pz_ep":^25}|'
        f'{"energy":^7}|'
        f'{"gen":^3}|'
        f'{"[simstat ]":^10}|'
        f'{"vertex x, y, z (mm)":^28}|'
        f'{"endpoint x, y, z (mm)":^28}|'
        f'{"mass":^7}|'
        f'{"charge":^6}|'
        f'{"[parents - daughters]"}'
    )
    print(get_sim_status_description())
    print(header)


def mc_print(mc):
    def mom_vec3(v):
        return f"{v[0]:^7.2f}, {v[1]:^7.2f}, {v[2]:^7.2f}"

    def pos_vec3(v):
        return f"{v[0]:^8.1f}, {v[1]:^8.1f}, {v[2]:^8.1f}"

    parents = ",".join(f"{p.id().index}" for p in mc.getParents())
    daughters = ",".join(f"{d.id().index}" for d in mc.getDaughters())

    try:
        name = Particle.from_pdgid(mc.getPDG()).name
    except:
        name = "unknown"
    mc_info = (
        f"{mc.id().index:^5d}|"
        f"{mc.getPDG():^10}|"
        f"{name:^15}|"
        f"{mom_vec3(mc.getMomentum())}|"
        f"{mom_vec3(mc.getMomentumAtEndpoint())}|"
        f"{mc.getEnergy():^7.3f}|"
        f"{mc.getGeneratorStatus():^3d}|"
        f"{get_sim_status(mc)}|"
        f"{pos_vec3(mc.getVertex())}|"
        f"{pos_vec3(mc.getEndpoint())}|"
        f"{mc.getMass():^7.3f}|"
        f"{mc.getCharge():^6.1f}|"
        f"[{parents}] - [{daughters}]"
    )
    print(mc_info)



def print_all_mcs(event):
    mc_print_header()
    for mc in event.get("MCParticles"):
        mc_print(mc)

if __name__ == "__main__":
    file = sys.argv[1]
    print_evt = int(sys.argv[2])
    reader = root_io.Reader(file)
    for idx, event in enumerate(reader.get("events")):
        if idx == print_evt:
            print_all_mcs(event)
            break