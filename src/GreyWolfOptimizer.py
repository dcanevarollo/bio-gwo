import random
from Blosum62 import Blosum62


class GreyWolfOptimizer:
    MATCH: int = 0
    GAP_B: int = 1
    GAP_A: int = 2
    GAP_OPEN: int = -5
    GAP_EXTEND: int = -1

    def __init__(self, population_size: int, dimension: int, max_iterations: int, sequence_a: str, sequence_b: str) -> None:
        self.population_size: int = population_size
        self.dimension: int = dimension
        self.max_iterations: int = max_iterations
        self.sequence_a: str = sequence_a
        self.sequence_b: str = sequence_b

        self.wolves: list[list[int]] = [[0] * dimension for _ in range(population_size)]
        for i in range(self.population_size):
            for j in range(self.dimension):
                self.wolves[i][j] = random.randint(0, 2)

        self.alpha: list[int] = [0] * dimension
        self.beta: list[int] = [0] * dimension
        self.delta: list[int] = [0] * dimension

        self.fitness: list[float] = [0.0] * population_size
        self.alpha_score: float = -1e9
        self.beta_score: float = -1e9
        self.delta_score: float = -1e9

        self.convergence: list[float] = []

    def optimize(self) -> None:
        for iteration in range(self.max_iterations):
            for i in range(self.population_size):
                self.fitness[i] = self.evaluate(self.wolves[i])

            self.update_hierarchy()
            self.update_positions(iteration)

            self.convergence.append(self.alpha_score)

    def print_best_alignment(self) -> None:
        best: list[int] = self.alpha

        aligned_a: list[str] = []
        aligned_b: list[str] = []

        i: int = 0
        j: int = 0
        for k in range(len(best)):
            operation: int = best[k]

            if operation == self.MATCH:
                if i < len(self.sequence_a) and j < len(self.sequence_b):
                    aligned_a.append(self.sequence_a[i])
                    aligned_b.append(self.sequence_b[j])
                    i += 1
                    j += 1
            elif operation == self.GAP_B:
                if i < len(self.sequence_a):
                    aligned_a.append(self.sequence_a[i])
                    aligned_b.append("-")
                    i += 1
            elif operation == self.GAP_A:
                if j < len(self.sequence_b):
                    aligned_a.append("-")
                    aligned_b.append(self.sequence_b[j])
                    j += 1

        print("\nBest alignment:")
        print("".join(aligned_a))
        print("".join(aligned_b))
        print("\nScore: " + str(self.alpha_score))

    def evaluate(self, wolf: list[int]) -> int:
        i: int = 0
        j: int = 0
        score: int = 0
        gap_a: bool = False
        gap_b: bool = False

        for k in range(len(wolf)):
            operation: int = wolf[k]

            if operation == self.MATCH:
                gap_a = False
                gap_b = False

                if i < len(self.sequence_a) and j < len(self.sequence_b):
                    score += Blosum62.score(self.sequence_a[i], self.sequence_b[j])
                    i += 1
                    j += 1
            elif operation == self.GAP_B:
                if not gap_b:
                    score += self.GAP_OPEN
                else:
                    score += self.GAP_EXTEND

                gap_b = True
                gap_a = False

                if i < len(self.sequence_a):
                    i += 1
            elif operation == self.GAP_A:
                if not gap_a:
                    score += self.GAP_OPEN
                else:
                    score += self.GAP_EXTEND

                gap_a = True
                gap_b = False

                if j < len(self.sequence_b):
                    j += 1

        return score

    def update_hierarchy(self) -> None:
        self.alpha_score = -1e9
        self.beta_score = -1e9
        self.delta_score = -1e9

        for i in range(self.population_size):
            if self.fitness[i] > self.alpha_score:
                self.delta_score = self.beta_score
                self.delta = self.beta.copy()

                self.beta_score = self.alpha_score
                self.beta = self.alpha.copy()

                self.alpha_score = self.fitness[i]
                self.alpha = self.wolves[i].copy()
            elif self.fitness[i] > self.beta_score:
                self.delta_score = self.beta_score
                self.delta = self.beta.copy()

                self.beta_score = self.fitness[i]
                self.beta = self.wolves[i].copy()
            elif self.fitness[i] > self.delta_score:
                self.delta_score = self.fitness[i]
                self.delta = self.wolves[i].copy()

    def update_positions(self, iteration: int) -> None:
        a: float = 2.0 - (2.0 * iteration / self.max_iterations)

        for i in range(self.population_size):
            for j in range(self.dimension):
                r1: float = random.random()
                r2: float = random.random()

                a1: float = 2 * a * r1 - a
                c1: float = 2 * r2
                d_alpha: float = abs(c1 * self.alpha[j] - self.wolves[i][j])
                x1: float = self.alpha[j] - a1 * d_alpha

                a2: float = 2 * a * random.random() - a
                c2: float = 2 * random.random()
                d_beta: float = abs(c2 * self.beta[j] - self.wolves[i][j])
                x2: float = self.beta[j] - a2 * d_beta

                a3: float = 2 * a * random.random() - a
                c3: float = 2 * random.random()
                d_delta: float = abs(c3 * self.delta[j] - self.wolves[i][j])
                x3: float = self.delta[j] - a3 * d_delta

                new_position: int = round((x1 + x2 + x3) / 3.0)

                self.wolves[i][j] = max(0, min(2, new_position))
