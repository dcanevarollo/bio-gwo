import random
import requests
import time
import tracemalloc
from io import StringIO
from gwo import GreyWolfOptimizer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def fetch_proteins(size=50) -> list[SeqRecord]:
    params = {
        "query": "(reviewed:true) AND (length:[100 TO 300]) AND (organism_id:9606)",
        "format": "fasta",
        "size": size
    }

    response = requests.get("https://rest.uniprot.org/uniprotkb/search", params=params)

    if response.status_code != 200:
        raise Exception("Error while accessing UnitProt")

    fasta_data = response.text
    sequences = list(SeqIO.parse(StringIO(fasta_data), "fasta"))

    return sequences

def sample_pairs(sequences: list[SeqRecord], n_pairs=10) -> list[tuple[SeqRecord, SeqRecord]]:
    pairs = []
    for _ in range(n_pairs):
        sequence1, sequence2 = random.sample(sequences, 2)
        pairs.append((sequence1, sequence2))

    return pairs

def run_alignment(seq1: str | SeqRecord, seq2: str | SeqRecord, population_size=100, max_iter=200) -> tuple[list[float], float, float, int]:
    if isinstance(seq1, SeqRecord):
        s1_name = seq1.name
        s1 = str(seq1.seq)
    else:
        s1_name = seq1
        s1 = seq1

    if isinstance(seq2, SeqRecord):
        s2_name = seq2.name
        s2 = str(seq2.seq)
    else:
        s2_name = seq2
        s2 = seq2
        
    print(f"Alinhando as seguintes proteínas:", s1_name, s2_name, sep="\n")

    optimizer = GreyWolfOptimizer(population_size, max_iter, s1, s2)

    tracemalloc.start()
    start_time = time.perf_counter()
    convergence, best_fitness = optimizer.run()
    end_time = time.perf_counter()
    _, memory_usage = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    elapsed_time = end_time - start_time

    optimizer.print_best_alignment()

    print(f"Tempo de execução: {elapsed_time:.2f} s")
    print(f"Memória consumida: {(memory_usage / (1024 * 1024)):.2f} MiB")

    return convergence, best_fitness, elapsed_time, memory_usage
