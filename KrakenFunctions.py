import sys
from ete3 import NCBITaxa
from collections import defaultdict

def parse_column4(col: str):
    """Parse k-mer assignations like 'taxid:count' into {taxid: count}, excluding 0 (unclassified)."""
    if not col.strip():
        return {}
    taxid_counts = {int(t.split(":")[0]): int(t.split(":")[1]) for t in col.split()}
    return {taxid: count for taxid, count in taxid_counts.items() if taxid != 0}

def get_lineage_dict(taxid,ncbi):
    """Return dict of {rank: name} for a given taxid."""
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)
    return {ranks[t]: names[t] for t in lineage if ranks[t] != "no rank"}

def get_LCA(taxids, ncbi):
    """Return the Lowest Common Ancestor (LCA) taxid of a list of taxids."""
    if not taxids:
        return None
    if len(taxids) == 1:
        return taxids[0]
    lineages = [set(ncbi.get_lineage(t)) for t in taxids]
    common = set.intersection(*lineages)
    if not common:
        return None
    # deepest common ancestor: iterate backwards over one lineage
    for taxid in ncbi.get_lineage(taxids[0])[::-1]:
        if taxid in common:
            return taxid
    return None

def classify_read(col: str, ncbi):
    """Classify a read from one side of a pair (k-mer histogram string)."""
    taxid_counts = parse_column4(col)
    if not taxid_counts:
        return {"LCA": None, "best_hit": None, "lineage": {}}

    # (a) LCA
    LCA_taxid = get_LCA(list(taxid_counts.keys()),ncbi)
    LCA_name = ncbi.get_taxid_translator([LCA_taxid]).get(LCA_taxid, "Unknown") if LCA_taxid else None

    # (b) Most supported genus/species
    best_taxid = max(taxid_counts, key=taxid_counts.get)
    lineage_dict = get_lineage_dict(best_taxid,ncbi)

    if "species" in lineage_dict:
        best_hit = (lineage_dict["species"], "species")
    elif "genus" in lineage_dict:
        best_hit = (lineage_dict["genus"], "genus")
    else:
        best_hit = (ncbi.get_taxid_translator([best_taxid]).get(best_taxid, "Unknown"), "other")

    return {
        "LCA": (LCA_taxid, LCA_name),
        "best_hit": best_hit,
        "lineage": lineage_dict
    }

def check_pair_consistency(r1, r2):
    """Check if paired reads are taxonomically coherent at the family level."""
    fam1 = r1["lineage"].get("family", None)
    fam2 = r2["lineage"].get("family", None)
    coherent = (fam1 is not None and fam1 == fam2)
    return coherent

def process_kraken2(input_file, output_file, ncbi):
    with open(input_file) as f, open(output_file, "w") as out:
        out.write("read_id\tLCA_R1\tbest_hit_R1\tfamily_R1\tLCA_R2\tbest_hit_R2\tfamily_R2\tpair_coherent\n")

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue

            read_id = parts[1]
            col5 = parts[4]
            if "A:" in col5:
                out.write(
                    f"{read_id}\t"
                    f"NA\t"
                    f"NA\t"
                    f"NA\t"
                    f"NA\t"
                    f"NA\t"
                    f"NA\t"
                    f"incoherent\n"
                )
                continue
            # split R1 / R2 by '|:|'
            elif "|:|" in col5:
                col_R1, col_R2 = col5.split("|:|")
            else:
                # if only one side provided, treat R2 as empty
                col_R1, col_R2 = col5, ""

            r1 = classify_read(col_R1.strip(),ncbi)
            r2 = classify_read(col_R2.strip(),ncbi)

            coherent = check_pair_consistency(r1, r2)

            out.write(
                f"{read_id}\t"
                f"{r1['LCA'][1] if r1['LCA'] else 'NA'}\t"
                f"{r1['best_hit'][0] if r1['best_hit'] else 'NA'}\t"
                f"{r1['lineage'].get('family','NA')}\t"
                f"{r2['LCA'][1] if r2['LCA'] else 'NA'}\t"
                f"{r2['best_hit'][0] if r2['best_hit'] else 'NA'}\t"
                f"{r2['lineage'].get('family','NA')}\t"
                f"{coherent}\n"
            )

def extract_taxid_taxline(name, ncbi):
    taxid = ncbi.get_name_translator([name])
    if name in taxid:
        tid = taxid[name][0]
        # Get lineage for the TaxID
        lineage = ncbi.get_lineage(tid)
        # Translate lineage TaxIDs back to names
        names = ncbi.get_taxid_translator(lineage)
        return taxid, lineage, names