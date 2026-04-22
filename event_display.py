#!/usr/bin/env python

from podio import root_io
import edm4hep
from kbhit import KBHit
import sys
import os
import argparse
import numpy as np
import cppyy
from itertools import cycle
import math

from dump_mc_table import mc_print_header, mc_print, print_all_mcs

colors = [
    "#D53E4F",  # red
    "#3288BD",  # blue
    "#FDAE61",  # orange
    "#5E4FA2",  # deep purple
    "#66C2A5",  # teal
    "#D01C8B",  # magenta
    "#F46D43",  # orange-red
    "#9970AB",  # medium purple
    "#018571",  # dark teal
    "#80CDC1",  # light teal
    "#9E0142",  # deep crimson
    "#4393C3",  # steel blue
    "#FFC107",  # amber
    "#762A83",  # dark violet
    "#7FBC41",  # lime green
    "#C51B7D",  # cerise
    "#313695",  # dark navy
    "#D6604D",  # salmon-red
    "#4D9221",  # olive green
    "#A50026",  # blood red
    "#006837",  # forest green
    "#2166AC",  # dark sky blue
    "#E08214",  # ochre
    "#35978F",  # slate teal
    "#E9A3C9",  # dusty rose
    "#BF812D",  # bronze
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
    parser.add_argument("--color-by", choices=["mc", "pandora", "labels"], default="mc",
                        help="Initial coloring mode: MC truth (default), Pandora PFOs, or numpy labels")
    parser.add_argument("--labels", metavar="FILE",
                        help="Path to .npz file with hit labels (index, collectionID, labels_hitpf)")
    parser.add_argument("--no-tracks", action="store_true",
                        help="Do not draw tracks, only hits")
    parser.add_argument("--min-hit-energy", type=float, default=0.0, metavar="GeV",
                        help="Skip hits with energy below this threshold (default: 0)")
    parser.add_argument("--show-mc", type=int, nargs="+", metavar="ID", default=[],
                        help="Only show hits/tracks from these MC truth IDs (e.g. --show-mc 23 42)")
    args = parser.parse_args()

    if args.color_by == "labels" and not args.labels:
        parser.error("--color-by labels requires --labels FILE")

    compact_path = args.compact_path
    file = args.file

    # Load numpy labels if provided.
    # Key: (rounded_collectionID, per_collection_index). collectionID is stored as float32
    # so we round through float32 on both sides to ensure keys match.
    hit2label = {}
    if args.labels:
        npz = np.load(args.labels)
        hit2label = {
            (int(cid), int(idx)): int(lbl)
            for cid, idx, lbl in zip(
                npz['collectionID'].astype(np.int64), # collection ID (float32-rounded)
                npz['index'].astype(np.int64),        # per-collection index
                npz['labels_hitpf'],
            )
        }

    available_modes = ["mc", "pandora"] + (["labels"] if args.labels else [])
    color_mode = args.color_by
    no_tracks = args.no_tracks
    min_hit_energy = args.min_hit_energy
    show_mc_ids = set(args.show_mc)

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

    # Override detector visualization colors.
    # drawDD4hepDetector reads colors from each DetElement's VisAttr, so setting
    # them here once means all subsequent draw calls use the custom colors.
    _dd4hep = cppyy.gbl.dd4hep
    def _set_det_color(det_element, r, g, b):
        try:
            vis = _dd4hep.VisAttr(str(det_element.name()) + "_custom_vis")
            vis.setColor(1.0, r / 255.0, g / 255.0, b / 255.0)
            det_element.setVisAttributes(vis)
        except Exception as e:
            print(f"Warning: could not set detector color for {det_element.name()}: {e}")

    for name_det_pair in detector.detectors():
        name = str(name_det_pair.first)
        nl   = name.lower()
        det_el = detector.detector(name)  # DetElement handle (has setVisAttributes)
        if "ecal" in nl:
            _set_det_color(det_el, 0xEE, 0xE7, 0xE2)   # #EEE7E2
        elif "hcal" in nl:
            _set_det_color(det_el, 26, 152, 80)          # rgb(26,152,80)
        elif any(t in nl for t in ["vertex", "vxd", "tracker", "sit", "ftd", "set", "tpc"]):
            _set_det_color(det_el, 255, 255, 255)         # white

    reader = root_io.Reader(file)
    kb = KBHit()

    def draw_event(event, evt_idx, color_mode):
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

        # Always build MC color map — needed for both display modes
        mcs = set()
        mcs.update(track2mc.values())
        mcs.update(hit2mc.values())
        mcs = sorted(mcs, key=lambda mc: mc.getEnergy(), reverse=True)
        mc2color, mc2id, id2mc = get_mcs_info(mcs)

        # Precompute allowed hit/track sets for --show-mc filter.
        # IDs refer to the MC particle's index in the MCParticles collection.
        if show_mc_ids:
            allowed_mcs = {mc for mc in mcs if mc.getObjectID().index in show_mc_ids}
            allowed_hits = {hit for hit, mc in hit2mc.items() if mc in allowed_mcs}
            allowed_tracks = {track for track, mc in track2mc.items() if mc in allowed_mcs}
        else:
            allowed_hits = allowed_tracks = None

        if color_mode == "pandora":
            pfos = list(event.get(PANDORA_PFO_COL))
            obj2id = {pfo: idx+1 for idx, pfo in enumerate(pfos)}
            id2obj = {idx+1: pfo for idx, pfo in enumerate(pfos)}

            # Color each PFO by the MC particle with the most hits inside it;
            # PFOs with no resolvable MC get a new unique color from the palette
            unresolved_colors = cycle(colors)
            obj2color = {}
            for pfo in pfos:
                hit_counts = {}
                for cluster in pfo.getClusters():
                    for hit in cluster.getHits():
                        if hit in hit2mc:
                            mc = hit2mc[hit]
                            hit_counts[mc] = hit_counts.get(mc, 0) + 1
                for track in pfo.getTracks():
                    for tracker_hit in track.getTrackerHits():
                        if tracker_hit in hit2mc:
                            mc = hit2mc[tracker_hit]
                            hit_counts[mc] = hit_counts.get(mc, 0) + 1
                if hit_counts:
                    dominant_mc = max(hit_counts, key=hit_counts.get)
                    obj2color[pfo] = mc2color.get(dominant_mc, next(unresolved_colors))
                else:
                    obj2color[pfo] = next(unresolved_colors)

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

        elif color_mode == "labels":
            # Map hits/tracks to their label via ObjectID lookup.
            # The .npz collectionIDs were saved as float32, so round through float32
            # to get the same value used as the key when building hit2label.
            def label_key(oid):
                return (int(np.float32(oid.collectionID)), oid.index)

            hit_assoc = {}
            for col in TRACKER_HIT_COLS + CALO_HIT_COLS:
                if col not in event.getAvailableCollections():
                    continue
                for hit in event.get(col):
                    key = label_key(hit.getObjectID())
                    if key in hit2label:
                        hit_assoc[hit] = hit2label[key]

            track_assoc = {}
            for track in event.get(TRACKS_COL):
                key = label_key(track.getObjectID())
                if key in hit2label:
                    track_assoc[track] = hit2label[key]

            # Color each label by its dominant MC (most hits), new color if unresolvable
            unique_labels = sorted(set(hit2label.values()))
            label_mc_counts = {}
            for hit, label in hit_assoc.items():
                if hit in hit2mc:
                    mc = hit2mc[hit]
                    label_mc_counts.setdefault(label, {})
                    label_mc_counts[label][mc] = label_mc_counts[label].get(mc, 0) + 1

            print(f"[labels] matched {len(hit_assoc)} hits and {len(track_assoc)} tracks out of {len(hit2label)} label entries")

            unresolved_colors = cycle(colors)
            label2color = {}
            for label in unique_labels:
                if label == 0:
                    label2color[label] = "#808080"
                elif label in label_mc_counts:
                    dominant_mc = max(label_mc_counts[label], key=label_mc_counts[label].get)
                    label2color[label] = mc2color.get(dominant_mc, next(unresolved_colors))
                else:
                    label2color[label] = next(unresolved_colors)

            obj2color = label2color
            obj2id = {label: idx+1 for idx, label in enumerate(unique_labels)}
            id2obj = {idx+1: label for idx, label in enumerate(unique_labels)}

        else:  # mc
            obj2color, obj2id, id2obj = mc2color, mc2id, id2mc
            track_assoc = track2mc
            hit_assoc = hit2mc

        if not no_tracks:
            for track in event.get(TRACKS_COL):
                if allowed_tracks is not None and track not in allowed_tracks:
                    continue
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
            is_tracker = col in TRACKER_HIT_COLS
            for hit in event.get(col):
                if allowed_hits is not None and hit not in allowed_hits:
                    continue
                energy = hit.getEDep() if is_tracker else hit.getEnergy()
                if energy < min_hit_energy:
                    continue
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

        color_mode_label = {"mc": "MC truth", "pandora": "Pandora PFO", "labels": "Labels"}[color_mode]
        print_header = f'Event # {evt_idx} [{color_mode_label} coloring]. Enter - next. ESC - quit. "c" - cycle coloring. Double click - print info. "p" - print all MCs.'
        print(f"{print_header:-^180}\n")
        if color_mode == "mc":
            mc_print_header()

        return id2obj

    for evt_idx, event in enumerate(reader.get("events")):
        if evt_idx==71:
            id2obj = draw_event(event, evt_idx, color_mode)

            while True:
                if kb.kbhit():
                    c = kb.getch()
                    if ord(c) == 27 or ord(c) == 10: # ESC or Enter has been hit
                        break
                    elif ord(c) == ord('c'):
                        color_mode = available_modes[(available_modes.index(color_mode) + 1) % len(available_modes)]
                        id2obj = draw_event(event, evt_idx, color_mode)
                    elif ord(c) == ord('p'):
                        print_all_mcs(event)

                selected_id = ced_selected_id_noblock()
                if selected_id >= 0 and selected_id in id2obj:
                    if color_mode == "pandora":
                        pfo = id2obj[selected_id]
                        print(f"Selected PFO #{selected_id}: E={pfo.getEnergy():.2f} GeV, PDG={pfo.getPDG()}, "
                            f"tracks={pfo.tracks_size()}, clusters={pfo.clusters_size()}")
                    elif color_mode == "labels":
                        print(f"Selected label #{id2obj[selected_id]}")
                    else:
                        mc_print(id2obj[selected_id])
                elif selected_id >= 0 and selected_id not in id2obj:
                    print("Selected hit has no assigned object!")
            kb.set_normal_term()

            if ord(c) == 27:
                exit()
