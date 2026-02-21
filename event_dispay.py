#!/usr/bin/env python

from podio import root_io
import edm4hep
from kbhit import KBHit
import sys
import os
import cppyy
from itertools import cycle
import math

from dump_mc_table import mc_print_header, mc_print, print_all_mcs

# created with Polychrome https://www.jstatsoft.org/article/view/v090c01
# R> set.seed(567629)
# R> createPalette(50, c("#2A95E8", "#E5629C"), range = c(5, 85), M = 100000)
colors = [
    "#2296E9", "#E3609A", "#0DF516", "#F5CB0D", "#FF0D00", "#0D4022", "#DB00FD", "#00EAE3",
    "#DFCBEC", "#F67A22", "#ADE391", "#6D0D89", "#4D16FE", "#773D0D", "#FE16D4", "#FC167E",
    "#F4CDA8", "#FE98FE", "#385F83", "#B182FE", "#562245", "#737100", "#9B1616", "#83B3A2",
    "#0DF1AE", "#8AD2F8", "#C7DE0D", "#F59284", "#007200", "#0D22B4", "#C48596", "#A30D88",
    "#D82A22", "#C0810D", "#32BE38", "#A28AD1", "#EA009B", "#807565", "#83163D", "#B20DBD",
    "#D5CD6C", "#2E7DFE", "#FF6882", "#94EB00", "#000D7E", "#B16AAA", "#9D22ED", "#FF6DD4",
    "#008E95", "#000D45"
]


def load_dependencies():
    for path in os.environ["CMAKE_PREFIX_PATH"].split(":"):
        inlcude_dir = os.path.join(path, "include")
        if os.path.isdir(inlcude_dir):
            cppyy.add_include_path( inlcude_dir )
        if "/marlinutil/" in path:
            cppyy.add_include_path( os.path.join(path, "include", "marlinutil") )

    cppyy.include("DDMarlinCED.h")

def init_event_display(host, port):
    ced_client_init(host, port)
    ced_register_elements()

def load_detector_geometry(compact_path):
    det_compact_file = os.path.join( compact_path )
    detector = cppyy.gbl.dd4hep.Detector.getInstance()
    detector.fromCompact(det_compact_file)
    return detector

def get_hit_mc_map(event, hit_rel_col):
    hit2mc = {}
    hit2weight = {}
    
    for link in event.get(hit_rel_col):
        if not link.getFrom().isAvailable() or not link.getTo().isAvailable():
            continue

        if isinstance(link.getTo(), edm4hep.MCParticle):
            mc = link.getTo()
        else:
            mc = link.getTo().getParticle()

        if link.getFrom() not in hit2mc:
            hit2mc[link.getFrom()] = mc
            hit2weight[link.getFrom()] = link.getWeight()
        else:
            if hit2weight[link.getFrom()] < link.getWeight():
                hit2mc[link.getFrom()] = mc
                hit2weight[link.getFrom()] = link.getWeight()
    return hit2mc

def get_mcs_info(mcs):
    mc2color = {mc : color for mc, color in zip(mcs, cycle(colors))}
    mc2id = {mc : id+1 for id, mc in enumerate(mcs)}
    id2mc = {id+1 : mc for id, mc in enumerate(mcs)}
    return mc2color, mc2id, id2mc

def get_track_momentum(track, bz):
    omega = track.omega
    if omega == 0.:
        return (0., 0., 0.)
    phi = track.phi
    tanL = track.tanLambda
    c_light = 299.792458 # mm/ns
    pt = (1e-6 * c_light * bz) / abs(omega)
    return (pt * math.cos(phi), pt * math.sin(phi), pt * tanL)

if __name__ == "__main__":
    compact_path = sys.argv[1]
    file = sys.argv[2]

    MC_PARTICLE_COL = "MCParticles"
    PANDORA_PFO_COL = "PandoraPFOs"
    CLUSTERS_COL = "PandoraClusters"
    CALOHIT_TO_MC_LINK_COL = "CalohitMCTruthLink"

    if "/ILD/" in compact_path:
        B_FIELD = 3.5
        TRACKS_COL = "MarlinTrkTracks"
        TRACK_TO_MC_LINK_COL = "MCTruthMarlinTrkTracksLink"
        TRACKER_HIT_COLS = ["VXDTrackerHits", "SITTrackerHits", "FTDPixelTrackerHits", "FTDStripTrackerHits", "SETTrackerHits", "TPCTrackerHits"]
        TRACKER_HIT_RELATION_COLS = ["VXDTrackerHitRelations", "FTDPixelTrackerHitRelations", "FTDStripTrackerHitRelations", "SITTrackerHitRelations", "SETTrackerHitRelations", "TPCTrackerHitRelations"]
        CALO_HIT_COLS = ["EcalBarrelCollectionRec", "EcalEndcapRingCollectionRec", "EcalEndcapsCollectionRec", "HcalBarrelCollectionRec", "HcalEndcapRingCollectionRec", "HcalEndcapsCollectionRec", "MUON", "LCAL", "LHCAL", "BCAL"]

    elif "/CLD/" in compact_path:
        B_FIELD = 2.0
        TRACKS_COL = "SiTracks"
        TRACK_TO_MC_LINK_COL = "SiTracksMCTruthLink"
        # TRACKER_HIT_COLS = ["VXDTrackerHits", "VXDEndcapTrackerHits", "ITrackerHits", "OTrackerHits"]
        # TRACKER_HIT_RELATION_COLS = ["VXDTrackerHitRelations", "VXDEndcapTrackerHitRelations"]
        TRACKER_HIT_COLS = []
        TRACKER_HIT_RELATION_COLS = []
        CALO_HIT_COLS = ["ECALBarrel", "ECALEndcap", "HCALBarrel", "HCALEndcap", "HCALOther", "MUON"]
    else:
        raise ValueError(f"Unknown detector in compact_path: {compact_path}")


    load_dependencies()
    from cppyy.gbl import ced_client_init, ced_register_elements, ced_new_event, ced_send_event, ced_hit_ID, CED_HIT_POINT, ced_selected_id_noblock
    std = cppyy.gbl.std
    DDMarlinCED = cppyy.gbl.DDMarlinCED


    ced_proc = init_event_display("localhost", 7286)
    detector = load_detector_geometry(compact_path)


    reader = root_io.Reader(file)
    kb = KBHit()
    for evt_idx, event in enumerate(reader.get("events")):
        ced_new_event()
        DDMarlinCED.drawDD4hepDetector( detector, False, std.vector[std.string]() )

        track2mc = {}
        for link in event.get(TRACK_TO_MC_LINK_COL):
            if not link.getFrom().isAvailable() or not link.getTo().isAvailable():
                continue
            track2mc[link.getFrom().id().index] = link.getTo()

        hit2mc = {}
        for hit_rel_col in TRACKER_HIT_RELATION_COLS + [CALOHIT_TO_MC_LINK_COL]:
            hit2mc.update( get_hit_mc_map(event, hit_rel_col) )

        mcs = set()
        mcs.update( track2mc.values() )
        mcs.update( hit2mc.values() )
        mcs = sorted(mcs, key=lambda mc: mc.getEnergy(), reverse=True)
        mc2color, mc2id, id2mc = get_mcs_info(mcs)

        for track in event.get(TRACKS_COL):
            if track.id().index in track2mc:
                mc = track2mc[track.id().index]
                color = mc2color[mc]
                id = mc2id[mc]
            else: 
                # print(f"WARNING: Found track with unassigned truth MC. Drawing with gray color. But this should never happen!")
                color = "#808080"
                id = 0

            ts = [ts for ts in track.getTrackStates() if ts.location == edm4hep.TrackState.AtIP][0]
            x, y, z = ts.referencePoint.x, ts.referencePoint.y, ts.referencePoint.z
            px,py,pz = get_track_momentum(ts, B_FIELD)
            charge = 1 if ts.omega > 0 else -1
            r_max = 2000
            z_max = 2350
            DDMarlinCED.drawHelix(B_FIELD, charge, x, y, z, px, py, pz, 1, 5, int(color.lstrip("#"), 16), 0.0, r_max, z_max, id)

        hit_size = 5
        for col in TRACKER_HIT_COLS:
            for hit in event.get(col):
                if hit in hit2mc:
                    mc = hit2mc[hit]
                    color = mc2color[mc]
                    id = mc2id[mc]
                else: 
                    # print(f"WARNING: Found tracker hit with unassigned truth MC. Drawing with gray color. But this should never happen!")
                    color = "#808080"
                    id = 0
                pos = hit.getPosition()
                ced_hit_ID(pos.x, pos.y, pos.z, CED_HIT_POINT, 1, hit_size, int(color.lstrip("#"), 16), id)

        for col in CALO_HIT_COLS:
            for hit in event.get(col):
                if hit in hit2mc:
                    mc = hit2mc[hit]
                    color = mc2color[mc]
                    id = mc2id[mc]
                else: 
                    # print(f"WARNING: Found hit with unassigned truth MC for method # {ced_layer_idx}. Drawing with gray color.")
                    color = "#808080"
                    id = 0

                pos = hit.getPosition()
                ced_hit_ID(pos.x, pos.y, pos.z, CED_HIT_POINT, 1, hit_size, int(color.lstrip("#"), 16), id)

        ced_send_event()

        print_header = f'Event # {evt_idx}. Enter - next event. ESC - quit. Double click - print MC info. "p" - print all MCs.'
        print(f"{print_header:-^180}\n")
        mc_print_header()
        while True:
            if kb.kbhit():
                c = kb.getch()
                if ord(c) == 27 or ord(c) == 10: # ESC or Enter has been hit
                    break
                elif ord(c) == ord('p'):
                    print_all_mcs(event)


            selected_id = ced_selected_id_noblock()
            if selected_id >= 0 and selected_id in id2mc:
                mc_print( id2mc[selected_id] )
            elif selected_id >= 0 and selected_id not in id2mc:
                print("Selected hit has no assigned MC!")
        kb.set_normal_term()

        if ord(c) == 27:
            exit()
