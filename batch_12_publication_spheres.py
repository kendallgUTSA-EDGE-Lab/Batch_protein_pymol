from pymol import cmd
import os
import glob
import csv

# Publication-ready batch PyMOL pipeline
# whole protein = spheres
# non-DUF = gray50
# DUF = tv_blue
# catalytic residues = tv_orange
# transparent background
# ray 2400 x 2400 PNG

BASE_DIR = r"C:\Users\green\Documents\EnzyMM_PyMOL"
PDB_DIR = os.path.join(BASE_DIR, "pdb_files")
OUTPUT_DIR = os.path.join(BASE_DIR, "pymol_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

PROTEINS = [
    {"gene":"CG14269","duf":"DUF4791","uniprot":"Q9W4S5","duf_start":25,"duf_end":187,"residues":[103,107,120]},
    {"gene":"CG2336","duf":"DUF1487","uniprot":"Q9VI33","duf_start":97,"duf_end":311,"residues":[150,186,187]},
    {"gene":"CG14537","duf":"DUF4791","uniprot":"Q9VLW3","duf_start":31,"duf_end":198,"residues":[176,177,180]},
    {"gene":"CG13539","duf":"DUF1487","uniprot":"Q9W1T5","duf_start":37,"duf_end":259,"residues":[51,53,88]},
    {"gene":"CG12516","duf":"DUF1487","uniprot":"Q9VAW6","duf_start":82,"duf_end":298,"residues":[102,103,106]},
    {"gene":"CG5623","duf":"DUF1487","uniprot":"Q8T485","duf_start":35,"duf_end":251,"residues":[183,201,233]},
    {"gene":"CG12637","duf":"DUF1487","uniprot":"Q9W2S6","duf_start":22,"duf_end":241,"residues":[80,82,83]},
    {"gene":"CG33723","duf":"DUF1091","uniprot":"Q4ABG8","duf_start":74,"duf_end":159,"residues":[110,114,117]},
    {"gene":"CG14013","duf":"DUF4464","uniprot":"Q9VMP6","duf_start":11,"duf_end":240,"residues":[130,180,198]},
    {"gene":"CG3528","duf":"DUF4464","uniprot":"Q9VQD0","duf_start":2,"duf_end":233,"residues":[184,198,203,205]},
    {"gene":"CG2854","duf":"DUF4461","uniprot":"Q9W4Z5","duf_start":181,"duf_end":487,"residues":[393,394,397]},
    {"gene":"CG14017","duf":"DUF4464","uniprot":"Q9VMP5","duf_start":24,"duf_end":237,"residues":[90,93,224]},
]

def find_structure(uniprot):
    hits = []
    for pattern in (
        os.path.join(PDB_DIR, f"*{uniprot}*.pdb"),
        os.path.join(PDB_DIR, f"*{uniprot}*.cif"),
    ):
        hits.extend(glob.glob(pattern))
    if not hits:
        return None
    for version in ("model_v6", "model_v5", "model_v4"):
        preferred = sorted(x for x in hits if version in os.path.basename(x))
        if preferred:
            return preferred[0]
    return sorted(hits)[0]

if not os.path.isdir(PDB_DIR):
    raise RuntimeError(f"PDB folder not found: {PDB_DIR}")

summary_rows = []

for item in PROTEINS:
    gene = item["gene"]
    duf = item["duf"]
    uniprot = item["uniprot"]
    duf_start = item["duf_start"]
    duf_end = item["duf_end"]
    residues = item["residues"]

    structure_file = find_structure(uniprot)
    if structure_file is None:
        summary_rows.append([gene, duf, uniprot, "", "PDB_NOT_FOUND"])
        continue

    cmd.reinitialize()
    cmd.load(structure_file, "protein")

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)
    cmd.set("depth_cue", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("specular", 0.20)
    cmd.set("shininess", 15)
    cmd.set("ambient", 0.45)
    cmd.set("direct", 0.55)
    cmd.set("reflect", 0.10)
    cmd.set("sphere_quality", 4)

    cmd.hide("everything", "protein")
    cmd.show("spheres", "protein")
    cmd.set("sphere_scale", 1.0, "protein")
    cmd.set("sphere_transparency", 0.0, "protein")

    cmd.color("gray50", "protein")

    cmd.select("DUF_region", f"protein and polymer.protein and resi {duf_start}-{duf_end}")
    cmd.color("tv_blue", "DUF_region")

    residue_string = "+".join(str(x) for x in residues)
    cmd.select("catalytic_residues", f"protein and polymer.protein and resi {residue_string}")
    cmd.color("tv_orange", "catalytic_residues")

    cmd.orient("protein")
    cmd.zoom("protein", buffer=8)

    prefix = f"{gene}_{duf}_{uniprot}"
    out_png = os.path.join(OUTPUT_DIR, prefix + "_publication.png")
    cmd.ray(2400, 2400)
    cmd.png(out_png)

    out_pse = os.path.join(OUTPUT_DIR, prefix + "_publication.pse")
    cmd.save(out_pse)

    summary_rows.append([gene, duf, uniprot, structure_file, "SUCCESS"])

summary_file = os.path.join(OUTPUT_DIR, "batch_summary.tsv")
with open(summary_file, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene","duf","uniprot","structure_file","status"])
    writer.writerows(summary_rows)
