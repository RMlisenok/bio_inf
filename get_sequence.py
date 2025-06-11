import requests

genes = ["ENSG00000225940", "ENSG00000226119", "ENSG00000226397"]

for gene in genes:
    try:
        server = "https://rest.ensembl.org"
        ext = f"/sequence/id/{gene}?type=genomic"

        r = requests.get(server + ext, headers={"Content-Type": "text/x-fasta"})

        if not r.ok:
            print(f"Ошибка при запросе гена {gene}")
            continue

        with open(f"{gene}.fasta", "w") as f:
            f.write(r.text)
        print(f"Последовательность гена {gene} сохранена.")
    except Exception as e:
        print(f"Ошибка при обработке гена {gene}: {str(e)}")