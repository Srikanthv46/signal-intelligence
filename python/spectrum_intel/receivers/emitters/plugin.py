"""Pipeline plugin: drone emitter association. Which receiver-defined signals are one
radio (a hopping link seen on several channels) and which radios form one drone
system (association.py: frame-clock lattice + oscillator ppm fingerprints).
Adds report["emitters"] and signal.decoded.fields.emitter_ids."""
from __future__ import annotations

from spectrum_intel.decode.base import Plugin


class EmitterAssociation(Plugin):
    name = "emitters"
    protocol = "drone emitter association"

    def postprocess(self, ctx, signals):
        from . import association as A
        from .families import DRONE_FAMILIES, _assoc_models
        cap = ctx.cap
        em, owner = [], {}
        for s in signals:
            for e in s.features.get("_emissions", []):
                if e.family in DRONE_FAMILIES:
                    owner[id(e)] = s
                    em.append(e)
        if not em:
            return {"emitters": []}
        em.sort(key=lambda e: (e.start_us, e.center_mhz))
        A.measure_fingerprints(cap.iq, cap.fs, cap.fc, em)
        emitters = A.associate(em, cap.fc, models=_assoc_models())
        out = []
        for k, m in enumerate(emitters):
            sigs = []
            for i in m.members:
                s = owner[id(em[i])]
                if s not in sigs:
                    sigs.append(s)
            eid = f"E{k}"
            for s in sigs:               # one signal can hold several radios' bursts (Herelink control + video)
                ids = s.decoded.setdefault("fields", {}).setdefault("emitter_ids", [])
                if eid not in ids:
                    ids.append(eid)
            d = m.to_dict()
            d["emitter_id"] = eid
            d["operator_meaning"] = DRONE_FAMILIES.get(m.family, m.family)
            d["n_bursts"] = len(d.pop("members"))
            out.append(d)
        return {"emitters": out}

    def finalize_report(self, rep):
        for em in rep.get("emitters", []):
            em["signal_ids"] = [s["id"] for s in rep["signals"]
                                if em["emitter_id"] in s["decoded"]["fields"].get("emitter_ids", [])]


PLUGINS = [EmitterAssociation()]
