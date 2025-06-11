import os
import gzip
import matplotlib.pyplot as plt

INPUT_DIR = "genomes"
OUTPUT_DIR = "genomes/processed"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def calculate_gc_content(file_path):
    total_len = 0
    gc_count = 0

    with gzip.open(file_path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            gc_count += line.count("G") + line.count("C")
            total_len += len(line)

    if total_len == 0:
        return 0
    return (gc_count / total_len) * 100

def process_all_genomes():
    gc_values = {}
    found = False

    for file in os.listdir(INPUT_DIR):
        if file.endswith(".fa.gz"):
            found = True
            species_name = file.replace(".dna.toplevel.fa.gz", "").replace(".", "_")
            input_path = os.path.join(INPUT_DIR, file)
            gc_content = calculate_gc_content(input_path)

            output_csv = os.path.join(OUTPUT_DIR, f"{species_name}_gc_content.csv")
            with open(output_csv, "w") as out:
                out.write("Species,GC_Content\n")
                out.write(f"{species_name},{gc_content:.2f}\n")

            print(f"✅ Обработан: {species_name} — GC: {gc_content:.2f}%")
            gc_values[species_name] = gc_content

    if not found:
        print("❌ Не найдено ни одного .fa.gz файла.")
    return gc_values

def plot_gc_content(gc_values):
    if not gc_values:
        print("⚠️ Нет данных для построения графика.")
        return

    species = list(gc_values.keys())
    values = list(gc_values.values())

    plt.figure(figsize=(12, 6))
    plt.bar(species, values, color="skyblue")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("GC Content (%)")
    plt.title("GC-состав геномов")
    plt.tight_layout()
    plot_path = os.path.join(OUTPUT_DIR, "gc_content_plot.png")
    plt.savefig(plot_path)
    print(f"📊 График сохранён в: {plot_path}")

if __name__ == "__main__":
    gc_data = process_all_genomes()
    plot_gc_content(gc_data)
    print("\n✅ Все файлы обработаны.")
