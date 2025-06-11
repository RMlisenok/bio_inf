import os
import requests
from bs4 import BeautifulSoup

# Параметры
ENSEMBL_RELEASE = '113'
GENOME_DIR = './genomes'
SPECIES_LIST = [
    "homo_sapiens",
    "pan_troglodytes",
    "mus_musculus",
    "rattus_norvegicus",
    "canis_lupus_familiaris",
    "bos_taurus",
    "gallus_gallus",
    "danio_rerio",
    "drosophila_melanogaster",
    "caenorhabditis_elegans"
    "saccharomyces_cerevisiae"
]

def download_genome(species, filename_contains="dna.toplevel.fa.gz"):
    base_url = f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/fasta/{species}/dna/"
    output_dir = os.path.join(GENOME_DIR, species)
    os.makedirs(output_dir, exist_ok=True)

    try:
        print(f"\n🔍 Обработка: {species}")
        r = requests.get(base_url)
        r.raise_for_status()

        soup = BeautifulSoup(r.text, "html.parser")
        for link in soup.find_all("a"):
            href = link.get("href")
            if href and filename_contains in href:
                full_url = base_url + href
                print(f"⬇️  Скачивание: {full_url}")
                file_response = requests.get(full_url, stream=True)
                file_response.raise_for_status()

                local_path = os.path.join(output_dir, href)
                with open(local_path, 'wb') as f:
                    for chunk in file_response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"✅ Сохранено в: {local_path}")
                return
        print(f"⚠️  Файл '{filename_contains}' не найден для {species}")
    except Exception as e:
        print(f"❌ Ошибка при загрузке {species}: {e}")

def main():
    print("🚀 Начало загрузки геномов Ensembl Release", ENSEMBL_RELEASE)
    for species in SPECIES_LIST:
        download_genome(species)
    print("🏁 Загрузка завершена.")

if __name__ == "__main__":
    main()
