import time
from GreyWolfOptimizer import GreyWolfOptimizer


def main() -> None:
    population_size: int = 30
    dimension: int = 20
    max_iterations: int = 100
    sequence_a: str = "HEAGAWGHEE"
    sequence_b: str = "PAWHEAE"

    gwo: GreyWolfOptimizer = GreyWolfOptimizer(population_size, dimension, max_iterations, sequence_a, sequence_b)

    start: float = time.time()
    gwo.optimize()
    end: float = time.time()

    print(f"Execution time: {(end - start) * 1000:.0f}ms")

    gwo.print_best_alignment()


if __name__ == "__main__":
    main()
