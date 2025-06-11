import sys
from Bio import SeqIO


def main():
    if len(sys.argv) != 4:
        print("Использование: python filter_extract.py <genome.fa> <input.bed> <output.fasta>")
        sys.exit(1)

    genome_file, bed_file, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

    try:
        # Загружаем геном
        genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

        with open(bed_file) as f, open(output_file, "w") as out:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                name, score, strand = parts[3], parts[4], parts[5]

                if chrom in genome:
                    seq = genome[chrom].seq[start:end]
                    if strand == '-':
                        seq = seq.reverse_complement()
                    out.write(f">{name}|{chrom}:{start}-{end}\n{seq}\n")
                else:
                    print(f"Предупреждение: {chrom} не найдена в геноме")

        print(f"Успешно обработан {bed_file}")

    except Exception as e:
        print(f"Ошибка: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()