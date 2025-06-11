import subprocess
from pathlib import Path
import shutil

def fix_bed_coordinates(bed_file: Path):
    fixed_lines = []
    with open(bed_file) as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if start > end:
                start, end = end, start
            parts[1], parts[2] = str(start), str(end)
            fixed_lines.append("\t".join(parts))

    # Создаём резервную копию
    backup = bed_file.with_suffix(".bed.bak")
    bed_file.rename(backup)
    with open(bed_file, "w") as out:
        out.write("\n".join(fixed_lines) + "\n")
    print(f"✅ Исправлен: {bed_file.name} (резерв: {backup.name})")

def run_bedtools_getfasta(bed_file: Path, genome_fasta: Path, output_fasta: Path):
    cmd = [
        "bedtools", "getfasta",
        "-fi", str(genome_fasta),
        "-bed", str(bed_file),
        "-s", "-name",
        "-fo", str(output_fasta)
    ]
    subprocess.run(cmd, check=True)

def count_fasta_sequences(fasta_path: Path) -> int:
    try:
        result = subprocess.run(
            ["grep", "-c", "^>", str(fasta_path)],
            capture_output=True, text=True, check=True
        )
        count = int(result.stdout.strip())
        print(f"📊 В FASTA {fasta_path.name}: {count} последовательностей")
        return count
    except subprocess.CalledProcessError:
        print(f"❌ Ошибка подсчёта последовательностей в {fasta_path.name}")
        return 0

def run_mafft_auto(input_fasta: Path, output_fasta: Path):
    print(f"🔗 MAFFT множественное выравнивание {input_fasta.name}...")
    try:
        with open(output_fasta, "w") as out_f:
            subprocess.run([
                "mafft", "--auto", str(input_fasta)
            ], stdout=out_f, check=True)
        print(f"✅ MAFFT завершён: {output_fasta.name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ MAFFT завершился с ошибкой: {e}")

def main():
    gene_id = "ENSG00000225940"
    results_dir = Path("nhmmer_results2")
    genome_dir = Path("genomes/processed1")
    final_out_dir = Path("final_outputs")
    final_out_dir.mkdir(exist_ok=True)

    input_hits_fa = final_out_dir / f"{gene_id}_hits.fa"
    aligned_fa = final_out_dir / f"{gene_id}_aligned.fa"

    genomes = [
        "homo_sapiens",
        "pan_troglodytes",
        "mus_musculus",
        "rattus_norvegicus",
        "canis_lupus_familiaris",
        "bos_taurus",
        "gallus_gallus",
        "danio_rerio"
    ]

    print(f"\n\U0001f9ec Обработка гена: {gene_id}")

    for species in genomes:
        print(f"📥 Обработка {species}...")
        bed_file = results_dir / f"{species.replace('_', ' ').title().replace(' ', '_')}.dna.toplevel.bed"
        genome_fasta = genome_dir / f"{species.replace('_', ' ').title().replace(' ', '_')}.dna.toplevel.fa"
        extracted_fa = results_dir / f"{species.replace('_', ' ').title().replace(' ', '_')}.dna.toplevel.extracted.fa"

        if not bed_file.exists():
            print(f"⚠️  Пропущено: {species} — нет .bed")
            continue

        try:
            fix_bed_coordinates(bed_file)
            run_bedtools_getfasta(bed_file, genome_fasta, extracted_fa)
        except subprocess.CalledProcessError as e:
            print(f"❌ Ошибка обработки {species}: {e}")

    if input_hits_fa.exists():
        count_fasta_sequences(input_hits_fa)
        run_mafft_auto(input_hits_fa, aligned_fa)
    else:
        print(f"⚠️  Файл {input_hits_fa} не найден")

if __name__ == "__main__":
    main()
