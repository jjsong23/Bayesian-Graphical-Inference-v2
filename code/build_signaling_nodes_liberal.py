#!/usr/bin/env python3
"""
LIBERAL / maximal-recall "all possible signaling participants" from a UniProt GAF + GO.

Difference from the precise build: broad enzyme roots (all kinases/phosphatases incl. lipid,
GTPases, cyclases, PDE, PI-PLC, NOS), receptors + ligands + regulators, adaptors/scaffolds,
second-messenger binding, AND the whole Biological-Process signaling subtree (any protein
annotated as participating in signal transduction / cell communication). Descendants are taken
over is_a AND part_of. No aspect filter — a protein qualifies via ANY matching MF or BP term.
Rows dropped only for NOT-negation and ND (no-data). All evidence kept; has_experimental flags it.

Ontology source: GO editors' file (go-edit.obo) — canonical go-basic is network-blocked here.
Asserted is_a/part_of only (no reasoner-inferred links), so sets are slightly conservative vs
official go-basic. Swap OBO_PATH for production.
"""
import obonet, networkx as nx, pandas as pd, os

OBO_PATH = "go.obo"
GAF_PATH = "/mnt/user-data/uploads/MOUSE-uniprot.gaf"
OUT_DIR  = "/mnt/user-data/outputs"

# functional class -> seed MF/BP root term(s); each expands to {root} + is_a/part_of descendants
ROOTS = {
 "receptor":                   ["GO:0060089","GO:0038023"],            # molecular transducer / signaling receptor activity
 "ligand":                     ["GO:0048018","GO:0005102"],            # receptor ligand activity / signaling receptor binding
 "receptor_regulator":         ["GO:0030545"],                          # signaling receptor regulator activity
 "kinase":                     ["GO:0016301"],                          # kinase activity (protein+lipid+small-mol)
 "phosphatase":                ["GO:0016791"],                          # phosphatase activity (protein+lipid)
 "gtpase":                     ["GO:0003924"],                          # GTPase activity (Ras/Rho/Rab/Arf switches)
 "gtpase_regulator":           ["GO:0030695","GO:0005085","GO:0005096"],# GTPase regulator / GEF / GAP
 "cyclase":                    ["GO:0004016","GO:0004383"],            # adenylate + guanylate cyclase
 "phosphodiesterase":          ["GO:0004112"],                          # cyclic-nucleotide PDE
 "phospholipase":              ["GO:0120548"],                          # phosphatidylinositol phospholipase C
 "nos":                        ["GO:0004517"],                          # nitric-oxide synthase
 "adaptor_scaffold":           ["GO:0060090","GO:0140378"],            # molecular adaptor / protein complex scaffold
 "kinase_phosphatase_binding": ["GO:0019900","GO:0019902","GO:0051018"],# kinase/phosphatase/PKA binding
 "second_messenger_binding":   ["GO:0030551","GO:0005516"],            # cyclic nucleotide binding / calmodulin binding
 "signaling_process":          ["GO:0023052","GO:0007165","GO:0035556","GO:0007154"], # BP signaling subtree
 "signaling_regulation":       ["GO:0009966"],                          # regulation of signal transduction
}
EXP_CODES = {"EXP","IDA","IPI","IMP","IGI","IEP","HTP","HDA","HMP","HGI","HEP"}

# ---- ontology ----------------------------------------------------------------------
g = obonet.read_obo(OBO_PATH)
sub = nx.MultiDiGraph(); sub.add_nodes_from(g.nodes(data=True))
sub.add_edges_from((u,v,k) for u,v,k in g.edges(keys=True) if k in ("is_a","part_of"))
name = lambda t: g.nodes.get(t,{}).get("name",t)
ns   = lambda t: g.nodes.get(t,{}).get("namespace","")
def alts(t):
    a=g.nodes[t].get("alt_id",[]); return [a] if isinstance(a,str) else list(a)

term2classes, seed_terms = {}, set()
for cls, roots in ROOTS.items():
    accepted=set()
    for r in roots:
        seed_terms.add(r); accepted.add(r); accepted |= nx.ancestors(sub, r)  # ancestors() = descendants in obonet
    for t in list(accepted):
        for tid in [t]+alts(t):
            term2classes.setdefault(tid,set()).add(cls)

# ---- GAF ---------------------------------------------------------------------------
cols=["DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID","DB_Reference","Evidence_Code",
      "With_From","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","Taxon","Date",
      "Assigned_By","Annotation_Extension","Gene_Product_Form_ID"]
gaf=pd.read_csv(GAF_PATH,sep="\t",comment="!",header=None,names=cols,dtype=str,keep_default_na=False)
n0=len(gaf)
gaf=gaf[~gaf.Qualifier.str.startswith("NOT")]
gaf=gaf[gaf.Evidence_Code!="ND"]
gaf=gaf[gaf.GO_ID.isin(term2classes)].copy()
gaf["classes"]=gaf.GO_ID.map(term2classes)

grp=gaf.groupby("DB_Object_ID")
nodes=pd.DataFrame({
 "symbol":           grp.DB_Object_Symbol.first(),
 "name":             grp.DB_Object_Name.first(),
 "classes":          grp.classes.agg(lambda c:";".join(sorted(set().union(*c)))),
 "n_classes":        grp.classes.agg(lambda c:len(set().union(*c))),
 "aspects":          grp.Aspect.agg(lambda s:";".join(sorted(set(s)))),
 "matched_go_ids":   grp.GO_ID.agg(lambda s:";".join(sorted(set(s)))),
 "n_matched_terms":  grp.GO_ID.agg(lambda s:len(set(s))),
 "evidence_codes":   grp.Evidence_Code.agg(lambda s:";".join(sorted(set(s)))),
 "has_experimental": grp.Evidence_Code.agg(lambda s:bool(set(s)&EXP_CODES)),
 "n_annotations":    grp.size(),
}).reset_index().rename(columns={"DB_Object_ID":"uniprot"}).sort_values("symbol").reset_index(drop=True)

os.makedirs(OUT_DIR,exist_ok=True)
nodes.to_csv(f"{OUT_DIR}/mouse_signaling_nodes_liberal.tsv",sep="\t",index=False)

# ---- GO-terms-used export (full expanded list, seed flagged) -----------------------
rows=[]
for tid,cs in term2classes.items():
    if tid in g.nodes:  # skip pure alt_ids in the readable export
        rows.append({"go_id":tid,"name":name(tid),"aspect":ns(tid),
                     "classes":";".join(sorted(cs)),"is_seed_root":tid in seed_terms})
terms=pd.DataFrame(rows).sort_values(["classes","aspect","go_id"]).reset_index(drop=True)
terms.to_csv(f"{OUT_DIR}/go_terms_used.tsv",sep="\t",index=False)

# ---- summary -----------------------------------------------------------------------
print(f"GAF rows {n0:,} -> matched (after NOT/ND drop) {len(gaf):,}")
print(f"\n=== LIBERAL NODE UNIVERSE: {len(nodes):,} gene products ===")
print(f"experimental support: {nodes.has_experimental.sum():,} ({100*nodes.has_experimental.mean():.0f}%)")
asp=nodes.aspects.value_counts()
print("qualified via aspect(s):", dict(asp))
print(f"seed root terms: {len(seed_terms)} | total distinct GO terms in net: {terms.go_id.nunique():,}")
print("\nper class (proteins; overlapping):")
for cls in ROOTS:
    m=nodes.classes.str.split(";").apply(lambda c:cls in c)
    print(f"  {cls:26s} {m.sum():5d}")
print("\nspot checks:")
for sym in ["Avpr2","Adcy6","Pde4d","Gnas","Prkaca","Akap5","Wnk1","Rgs2","Plcb1","Nos1","Camk2d","Arrb1"]:
    r=nodes[nodes.symbol.str.lower()==sym.lower()]
    if len(r):
        row=r.iloc[0]; print(f"  {sym:8s} {row.uniprot:11s} [{row.classes}]")
    else:
        print(f"  {sym:8s} — not captured")
