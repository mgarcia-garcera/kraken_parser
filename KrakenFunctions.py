import sys
from ete3 import NCBITaxa
from collections import defaultdict

def parse_column4(col4: str):
    """Parse Kraken2 column 4 into {taxid: count}, excluding 0 (unclassified)."""
    taxid_counts = {int(t.split(":")[0]): int(t.split(":")[1]) for t in col4.split()}
    return {taxid: count for taxid, count in taxid_counts.items() if taxid != 0}

def get_lineage_dict(taxid, ncbi):
    """Return dict of {rank: name} for a given taxid."""
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    ranks = ncbi.get_rank(lineage)
    return {ranks[t]: names[t] for t in lineage if ranks[t] != "no rank"}

def get_LCA(taxids, ncbi):
    """Return the Lowest Common Ancestor (LCA) taxid of a list of taxids."""
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

def classify_read(col4, ncbi):
    """Classify a read given Kraken2 col4 string."""
    taxid_counts = parse_column4(col4)
    if not taxid_counts:
        return {"LCA": None, "best_hit": None, "lineage": {}}

    # (a) LCA
    LCA_taxid = get_LCA(list(taxid_counts.keys()))
    LCA_name = ncbi.get_taxid_translator([LCA_taxid]).get(LCA_taxid, "Unknown") if LCA_taxid else None

    # (b) Most supported genus/species
    best_taxid = max(taxid_counts, key=taxid_counts.get)
    lineage_dict = get_lineage_dict(best_taxid,ncbi)

    best_hit = None
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

def process_kraken2(input_file, output_file,ncbi):
    pairs = defaultdict(dict)

    # Read Kraken2 output
    with open(input_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            read_id, classification, length, col4 = parts[0], parts[1], parts[2], parts[3]
            # Group by read ID prefix (before /1 or /2 if present)
            prefix = read_id.split()[0].replace("/1", "").replace("/2", "")
            if read_id.endswith("/1"):
                pairs[prefix]["R1"] = classify_read(col4,ncbi)
            elif read_id.endswith("/2"):
                pairs[prefix]["R2"] = classify_read(col4,ncbi)
            else:
                # single-end read
                pairs[prefix]["R1"] = classify_read(col4,ncbi)

    # Write TSV
    with open(output_file, "w") as out:
        out.write("read_id\tLCA\tbest_hit\tfamily\tpair_coherent\n")
        for rid, data in pairs.items():
            r1 = data.get("R1", None)
            r2 = data.get("R2", None)

            if r1 and r2:  # paired
                coherent = check_pair_consistency(r1, r2)
                fam = r1["lineage"].get("family", None) or r2["lineage"].get("family", None)
                out.write(f"{rid}\t{r1['LCA'][1] if r1['LCA'] else 'NA'}"
                          f"\t{r1['best_hit'][0] if r1['best_hit'] else 'NA'}"
                          f"\t{fam if fam else 'NA'}"
                          f"\t{coherent}\n")
            else:  # single-end
                r = r1 or r2
                fam = r["lineage"].get("family", None)
                out.write(f"{rid}\t{r['LCA'][1] if r['LCA'] else 'NA'}"
                          f"\t{r['best_hit'][0] if r['best_hit'] else 'NA'}"
                          f"\t{fam if fam else 'NA'}"
                          f"\tNA\n")