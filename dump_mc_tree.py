#! /usr/bin/env python3

from podio import root_io
from particle import Particle
from particle.particle.particle import ParticleNotFound
import sys
import graphviz
import math
from dump_mc_table import get_sim_status


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
    label += f"Gen/Sim: {mc.getGeneratorStatus()}/{get_sim_status(mc)}"
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

    dot.render(directory='output')


if __name__ == "__main__":
    file = sys.argv[1]
    print_evt = int(sys.argv[2])
    reader = root_io.Reader(file)
    for idx, event in enumerate(reader.get("events")):
        if idx == print_evt:
            draw_mc_tree(event, idx)
            break