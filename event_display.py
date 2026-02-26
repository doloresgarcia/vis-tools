#!/usr/bin/env python

from podio import root_io
import edm4hep
from kbhit import KBHit
import sys
import os
import argparse
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
    if hit_rel_col not in event.getAvailableCollections():
        print(f"Warning: Hit relation collection {hit_rel_col} not found in event.")
        return hit2mc

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
    parser = argparse.ArgumentParser(description="Event display for ILD/CLD detectors")
    parser.add_argument("compact_path", help="Path to detector compact file")
    parser.add_argument("file", help="Path to input ROOT file")
    parser.add_argument("--color-by", choices=["mc", "pandora"], default="mc",
                        help="Color hits and tracks by MC truth (default) or Pandora PFOs")
    args = parser.parse_args()

    compact_path = args.compact_path
    file = args.file
    use_pandora_coloring = args.color_by == "pandora"

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
        TRACKS_COL = "SiTracks_Refitted"
        TRACK_TO_MC_LINK_COL = "SiTracksMCTruthLink"
        TRACKER_HIT_COLS = ["VXDTrackerHits", "VXDEndcapTrackerHits", "ITrackerHits", "OTrackerHits"]
        TRACKER_HIT_RELATION_COLS = ["VXDTrackerHitRelations", "VXDEndcapTrackerHitRelations"]
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

    def draw_event(event, evt_idx, use_pandora):
        ced_new_event()
        subdet_names = std.vector[std.string]()
        for name_det_pair in detector.detectors():
            name = str(name_det_pair.first)
            if "yoke" not in name.lower() and "solenoid" not in name.lower():
                subdet_names.push_back(name)
        DDMarlinCED.drawDD4hepDetector(detector, False, subdet_names)

        track2mc = {}
        for link in event.get(TRACK_TO_MC_LINK_COL):
            if not link.getFrom().isAvailable() or not link.getTo().isAvailable():
                continue
            track2mc[link.getFrom()] = link.getTo()

        hit2mc = {}
        for hit_rel_col in TRACKER_HIT_RELATION_COLS + [CALOHIT_TO_MC_LINK_COL]:
            hit2mc.update(get_hit_mc_map(event, hit_rel_col))

        if use_pandora:
            pfos = list(event.get(PANDORA_PFO_COL))
            obj2color = {pfo: color for pfo, color in zip(pfos, cycle(colors))}
            obj2id = {pfo: idx+1 for idx, pfo in enumerate(pfos)}
            id2obj = {idx+1: pfo for idx, pfo in enumerate(pfos)}

            track_assoc = {}
            for pfo in pfos:
                for track in pfo.getTracks():
                    track_assoc[track] = pfo

            calohit2pfo = {}
            for pfo in pfos:
                for cluster in pfo.getClusters():
                    for hit in cluster.getHits():
                        calohit2pfo[hit] = pfo

            mc2track = {mc: track for track, mc in track2mc.items()}
            hit_assoc = dict(calohit2pfo)
            for hit, mc in hit2mc.items():
                if hit not in hit_assoc and mc in mc2track and mc2track[mc] in track_assoc:
                    hit_assoc[hit] = track_assoc[mc2track[mc]]
        else:
            mcs = set()
            mcs.update(track2mc.values())
            mcs.update(hit2mc.values())
            mcs = sorted(mcs, key=lambda mc: mc.getEnergy(), reverse=True)
            obj2color, obj2id, id2obj = get_mcs_info(mcs)
            track_assoc = track2mc
            hit_assoc = hit2mc

        for track in event.get(TRACKS_COL):
            if track in track_assoc:
                obj = track_assoc[track]
                color = obj2color[obj]
                id = obj2id[obj]
            else:
                color = "#808080"
                id = 0
            ts = [ts for ts in track.getTrackStates() if ts.location == edm4hep.TrackState.AtIP][0]
            x, y, z = ts.referencePoint.x, ts.referencePoint.y, ts.referencePoint.z
            px, py, pz = get_track_momentum(ts, B_FIELD)
            charge = 1 if ts.omega > 0 else -1
            DDMarlinCED.drawHelix(B_FIELD, charge, x, y, z, px, py, pz, 1, 5, int(color.lstrip("#"), 16), 0.0, 2000, 2350, id)

        hit_size = 5
        for col in TRACKER_HIT_COLS + CALO_HIT_COLS:
            if col not in event.getAvailableCollections():
                print(f"WARNING: hit collection {col} not found in event. Skipping.")
                continue
            for hit in event.get(col):
                if hit in hit_assoc:
                    obj = hit_assoc[hit]
                    color = obj2color[obj]
                    id = obj2id[obj]
                else:
                    color = "#808080"
                    id = 0
                pos = hit.getPosition()
                ced_hit_ID(pos.x, pos.y, pos.z, CED_HIT_POINT, 1, hit_size, int(color.lstrip("#"), 16), id)

        ced_send_event()

        color_mode_label = "Pandora PFO" if use_pandora else "MC truth"
        print_header = f'Event # {evt_idx} [{color_mode_label} coloring]. Enter - next. ESC - quit. "c" - toggle coloring. Double click - print info. "p" - print all MCs.'
        print(f"{print_header:-^180}\n")
        if not use_pandora:
            mc_print_header()

        return id2obj

    for evt_idx, event in enumerate(reader.get("events")):
        if evt_idx==71:
            id2obj = draw_event(event, evt_idx, use_pandora_coloring)

            while True:
                if kb.kbhit():
                    c = kb.getch()
                    if ord(c) == 27 or ord(c) == 10: # ESC or Enter has been hit
                        break
                    elif ord(c) == ord('c'):
                        use_pandora_coloring = not use_pandora_coloring
                        id2obj = draw_event(event, evt_idx, use_pandora_coloring)
                    elif ord(c) == ord('p'):
                        print_all_mcs(event)

                selected_id = ced_selected_id_noblock()
                if selected_id >= 0 and selected_id in id2obj:
                    if use_pandora_coloring:
                        pfo = id2obj[selected_id]
                        print(f"Selected PFO #{selected_id}: E={pfo.getEnergy():.2f} GeV, PDG={pfo.getPDG()}, "
                              f"tracks={pfo.tracks_size()}, clusters={pfo.clusters_size()}")
                    else:
                        mc_print(id2obj[selected_id])
                elif selected_id >= 0 and selected_id not in id2obj:
                    print("Selected hit has no assigned object!")
            kb.set_normal_term()

            if ord(c) == 27:
                exit()
