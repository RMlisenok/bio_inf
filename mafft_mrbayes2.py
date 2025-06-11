import os
import subprocess
from pathlib import Path
import shutil

GENES = {
    "ENSG00000226119": "nhmmer_results",
    "ENSG00000226397": "nhmmer_results3"
}

GENOMES = [
    "homo_sapiens",
    "pan_troglodytes",
    "mus_musculus",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "bos_taurus",
    "gallus_gallus",
    "danio_rerio"
]

PROJECT_DIR = Path("/home/lisa/projects/python_projects/genes")
QUERY_DIR = PROJECT_DIR / "query_sequences"
OUTPUT_DIR = PROJECT_DIR / "final_outputs_flexible"
OUTPUT_DIR.mkdir(exist_ok=True)

def fix_bed_coordinates(bed_file: Path):
    fixed_lines = []
    with open(bed_file) as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if start > end:
                start, end = end, start
            parts[1], parts[2] = str(start), str(end)
            fixed_lines.append("\t".join(parts))

    backup = bed_file.with_suffix(".bed.bak")
    bed_file.rename(backup)
    with open(bed_file, "w") as out:
        out.write("\n".join(fixed_lines) + "\n")
    print(f"✅ Исправлен: {bed_file.name} (резерв: {backup.name})")

def find_matching_file(folder: Path, genome: str, suffix: str):
    genome_lower = genome.lower()
    for file in folder.glob(f"*{suffix}"):
        if genome_lower in file.name.lower():
            return file
    return None

def clean_headers(input_fa: Path, output_fa: Path):
    with open(output_fa, "w") as out:
        for line in open(input_fa):
            if line.startswith(">"):
                clean = ''.join(c for c in line if c.isalnum() or c in ">\n_")
                out.write(clean)
            else:
                out.write(line)

def concatenate_fasta(fasta_list, output_path):
    with open(output_path, "w") as out:
        for fa in fasta_list:
            if fa.exists():
                out.write(open(fa).read())

def run_mafft(input_fa: Path, output_fa: Path):
    print("🔗 MAFFT множественное выравнивание...")
    subprocess.run(["mafft", "--auto", str(input_fa)], stdout=open(output_fa, "w"), stderr=subprocess.STDOUT)
    print(f"✅ Выравнивание сохранено: {output_fa.name}")

def process_gene(gene_id, nhmmer_folder):
    print(f"\n🧬 Обработка гена: {gene_id}")
    result_fasta = OUTPUT_DIR / f"{gene_id}_hits.fa"
    result_cleaned = OUTPUT_DIR / f"{gene_id}_hits_cleaned.fa"
    aligned = OUTPUT_DIR / f"{gene_id}_aligned.fa"
    query_fa = QUERY_DIR / f"{gene_id}.fasta"

    fasta_paths = []

    for genome in GENOMES:
        print(f"📥 Обработка {genome}...")
        bed = find_matching_file(PROJECT_DIR / nhmmer_folder, genome, ".bed")
        fa = find_matching_file(PROJECT_DIR / nhmmer_folder, genome, ".fa")

        if bed and fa:
            fix_bed_coordinates(bed)
            extracted = bed.with_suffix(".extracted.fa")
            subprocess.run([
                "bedtools", "getfasta",
                "-fi", str(fa),
                "-bed", str(bed),
                "-s", "-name",
                "-fo", str(extracted)
            ])
            fasta_paths.append(extracted)
        else:
            print(f"⚠️  Пропущено: {genome} — нет .bed или .fa")

    # Добавляем человека
    if query_fa.exists():
        fasta_paths.append(query_fa)
    else:
        print(f"⚠️  Нет исходной последовательности: {query_fa.name}")

    if not fasta_paths:
        print("❌ Нет доступных последовательностей")
        return

    concatenate_fasta(fasta_paths, result_fasta)
    clean_headers(result_fasta, result_cleaned)

    # Подсчёт последовательностей
    n_seqs = sum(1 for line in open(result_cleaned) if line.startswith(">"))
    print(f"📊 Всего последовательностей: {n_seqs}")

    if n_seqs >= 2:
        run_mafft(result_cleaned, aligned)
    else:
        print("⚠️  Недостаточно последовательностей для выравнивания (нужно минимум 2)")

def main():
    for gene, nhmmer_folder in GENES.items():
        process_gene(gene, nhmmer_folder)

if __name__ == "__main__":
    main()
