#! /usr/bin/env python3

from podio import root_io
from particle import Particle
from particle.particle.particle import ParticleNotFound
import sys
import graphviz
import math

def get_simulator_status_string(mc):
    if mc is None:
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


def get_pdg_name(pdg):
    name = "unknown"
    try:
        name = Particle.from_pdgid(pdg).name
    except ParticleNotFound:
        name = f"PDG {pdg}"
    return name

def get_mc_label(mc):
    label = f"{get_pdg_name(mc.getPDG())} / id {mc.id().index}\n"
    label += f"Vtx: ({mc.getVertex()[0]:.1f},{mc.getVertex()[1]:.1f},{mc.getVertex()[2]:.1f}) mm\n"
    label += f"End: ({mc.getEndpoint()[0]:.1f},{mc.getEndpoint()[1]:.1f},{mc.getEndpoint()[2]:.1f}) mm\n"
    label += f"Mom: ({mc.getMomentum()[0]:.1f},{mc.getMomentum()[1]:.1f},{mc.getMomentum()[2]:.1f}) GeV\n"
    label += f"E: {mc.getEnergy():.2f} GeV\n"
    label += f"Gen/Sim: {mc.getGeneratorStatus()}/{get_simulator_status_string(mc)}"
    return label

def draw_mc_tree(event, idx):
    mc2idx = { mc: idx for idx, mc in enumerate(event.get("MCParticles")) }

    dot = graphviz.Digraph(f"Event{idx}", comment='MC Particle Tree', format='svg', node_attr={'shape': 'ellipse','style': 'filled'})

    for idx, mc in enumerate(event.get("MCParticles")):
        if mc.getGeneratorStatus() > 1:
            col = 'dimgrey'
        elif mc.getGeneratorStatus() == 1:
            col = 'white'
        elif mc.getGeneratorStatus() == 0:
            col = 'lightgoldenrod1'


        dot.node(str(idx), get_mc_label(mc), fillcolor=col)
        for daughter in mc.getDaughters():
            dot.edge(str(idx), str(mc2idx[daughter]))

    dot.render(directory='doctest-output')


if __name__ == "__main__":
    file = sys.argv[1]
    print_evt = int(sys.argv[2])
    reader = root_io.Reader(file)
    for idx, event in enumerate(reader.get("events")):
        if idx == print_evt:
            draw_mc_tree(event, idx)
            break