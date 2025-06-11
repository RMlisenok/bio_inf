import json
from pathlib import Path
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO

def run_mafft(input_fasta: str, output_fasta: str):
    """Run MAFFT alignment and save result"""
    print(f"🔗 MAFFT выравнивание: {input_fasta}")
    mafft_cline = MafftCommandline(input=input_fasta)
    stdout, _ = mafft_cline()
    with open(output_fasta, "w") as f:
        f.write(stdout)
    print(f"✅ Сохранено выравнивание: {output_fasta}")

def tree_to_jvp(tree, output_path):
    """Convert Biopython tree to Jawline .jvp format"""
    def traverse(clade):
        node = {
            "name": str(clade.name) if clade.name else "",
            "length": clade.branch_length if clade.branch_length else 0,
            "children": []
        }
        for sub in clade.clades:
            node["children"].append(traverse(sub))
        return node

    data = {
        "format": "jvp",
        "tree": traverse(tree.root)
    }

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"📁 Дерево сохранено как .jvp: {output_path}")

def build_jvp_tree(input_fasta: str):
    input_path = Path(input_fasta)
    aligned_path = input_path.with_name(input_path.stem + "_aligned_tmp.fa")
    jvp_path = input_path.with_suffix(".jvp")

    if not input_path.exists():
        print(f"❌ Файл не найден: {input_path}")
        return

    # Step 1: Align
    run_mafft(str(input_path), str(aligned_path))

    # Step 2: Build tree
    alignment = AlignIO.read(str(aligned_path), "fasta")
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Step 3: Save in JVP format
    tree_to_jvp(tree, jvp_path)

# 🚀 Построение дерева
build_jvp_tree("ENSG00000225940_hits.fa")
