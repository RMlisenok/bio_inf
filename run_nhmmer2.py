import os
import subprocess
from pathlib import Path

# Путь к геномам и последовательности-запросу
GENOMES_DIR = Path("./genomes/processed1")
QUERY = Path("/home/lisa/projects/python_projects/genes/query_sequences/ENSG00000225940.fasta")
OUTPUT_DIR = Path("./nhmmer_results2")
OUTPUT_DIR.mkdir(exist_ok=True)

# E-value threshold
E_THRESHOLD = 1e-5

def run_nhmmer(genome_fasta, query_fasta):
    output_tsv = OUTPUT_DIR / (genome_fasta.stem + ".tsv")
    cmd = [
        "nhmmer",
        "--tblout", str(output_tsv),
        str(query_fasta),
        str(genome_fasta)
    ]
    print(f"🚀 Запуск nhmmer для {genome_fasta.name}")
    subprocess.run(cmd, check=True)
    return output_tsv

def parse_and_filter_tsv(tsv_path, e_threshold=1e-10):
    matches = []
    with open(tsv_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            try:
                target_name = fields[0]
                alifrom = int(fields[6])
                alito = int(fields[7])
                strand = fields[11]
                evalue = float(fields[12])
                print(f"Проверка: найдено совпадение с e-value={evalue} на {target_name}:{alifrom}-{alito} ({strand})")
            except Exception as ex:
                print(f"⚠️ Ошибка парсинга строки: {line.strip()}\n{ex}")
                continue

            if evalue <= e_threshold:
                matches.append((target_name, alifrom, alito, strand))
    return matches



def write_bed_file(matches, bed_path):
    """
    Записывает результаты в BED файл.
    BED координаты — 0-ориентированные, полузакрытые:
    start = alifrom - 1, end = alito
    """
    with open(bed_path, "w") as f:
        for chrom, start, end, strand in matches:
            bed_start = start - 1  # BED формат 0-based
            bed_end = end          # BED end - не включительно
            # Формат: chrom, start, end, name, score, strand
            # name и score можно оставить пустыми или задать фиктивные
            f.write(f"{chrom}\t{bed_start}\t{bed_end}\t.\t0\t{strand}\n")

def extract_sequences(genome_fasta, bed_path, fasta_out):
    """
    Извлекает последовательности с помощью bedtools getfasta.
    """
    cmd = [
        "bedtools", "getfasta",
        "-fi", str(genome_fasta),
        "-bed", str(bed_path),
        "-s",  # учитываем strand
        "-name",  # сохраняем имя из bed (хромосому)
        "-fo", str(fasta_out)
    ]
    print(f"🔍 Извлечение последовательностей из {genome_fasta.name}")
    subprocess.run(cmd, check=True)

def main():
    all_results = []
    for genome_fasta in GENOMES_DIR.glob("*.fa"):
        try:
            tsv_path = run_nhmmer(genome_fasta, QUERY)
            matches = parse_and_filter_tsv(tsv_path, e_threshold=E_THRESHOLD)
            print(f"✅ Найдено {len(matches)} совпадений в {genome_fasta.name}")

            if matches:
                bed_path = OUTPUT_DIR / f"{genome_fasta.stem}.bed"
                fasta_out = OUTPUT_DIR / f"{genome_fasta.stem}_extracted.fa"

                write_bed_file(matches, bed_path)
                extract_sequences(genome_fasta, bed_path, fasta_out)

                print(f"💾 BED файл сохранён: {bed_path}")
                print(f"💾 Последовательности для MAFFT сохранены: {fasta_out}")

            all_results.extend(matches)
        except Exception as e:
            print(f"❌ Ошибка с файлом {genome_fasta.name}: {e}")

    if not all_results:
        print("⚠️  Ничего не найдено по заданному e-value.")
    else:
        print(f"\n🎯 Всего подходящих совпадений: {len(all_results)}")

if __name__ == "__main__":
    main()
