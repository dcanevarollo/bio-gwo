import random
from Bio.Align import substitution_matrices


MATCH: int = 0
GAP_B: int = 1
GAP_A: int = 2
GAP_OPEN: int = -5
GAP_EXTEND: int = -1

class GreyWolfOptimizer:

    def __init__(self, population_size: int, max_iter: int, seq1: str, seq2: str) -> None:
        self.__population_size = population_size
        self.__dim = len(seq1) + len(seq2)
        self.__max_iter = max_iter
        self.__seq1 = seq1.upper()
        self.__seq2 = seq2.upper()

        self.__wolves = [[0] * self.__dim for _ in range(population_size)]

        random.seed(42)
        for i in range(self.__population_size):
            for j in range(self.__dim):
                # 70% MATCH, 15% GAP_B, 15% GAP_A for better initialization
                rand = random.random()
                if rand < 0.7:
                    self.__wolves[i][j] = 0  # MATCH
                elif rand < 0.85:
                    self.__wolves[i][j] = 1  # GAP_B
                else:
                    self.__wolves[i][j] = 2  # GAP_A

        self.__alpha = [0] * self.__dim
        self.__beta = [0] * self.__dim
        self.__delta = [0] * self.__dim

        self.__fitness = [0.0] * population_size
        self.__alpha_score = -1e9
        self.__beta_score = -1e9
        self.__delta_score = -1e9

        self.__blosum62 = substitution_matrices.load("BLOSUM62")

    def run(self) -> tuple[list[float], float]:
        convergence: list[float] = []

        for iteration in range(self.__max_iter):
            for i in range(self.__population_size):
                self.__fitness[i] = self.__evaluate(self.__wolves[i])

            self.__update_hierarchy()
            self.__update_positions(iteration)

            convergence.append(self.__alpha_score)

        return convergence, self.__alpha_score

    def print_best_alignment(self) -> None:
        best: list[int] = self.__alpha

        aligned_a: list[str] = []
        aligned_b: list[str] = []

        i: int = 0
        j: int = 0
        for k in range(len(best)):
            operation: int = best[k]

            if operation == MATCH:
                if i < len(self.__seq1) and j < len(self.__seq2):
                    aligned_a.append(self.__seq1[i])
                    aligned_b.append(self.__seq2[j])
                    i += 1
                    j += 1
            elif operation == GAP_B:
                if i < len(self.__seq1):
                    aligned_a.append(self.__seq1[i])
                    aligned_b.append("-")
                    i += 1
            elif operation == GAP_A:
                if j < len(self.__seq2):
                    aligned_a.append("-")
                    aligned_b.append(self.__seq2[j])
                    j += 1

        print("\nBest alignment:")
        print("".join(aligned_a))
        print("".join(aligned_b))
        print("\nScore: " + str(self.__alpha_score))

    def __evaluate(self, wolf: list[int]) -> float:
        i, j = 0, 0
        score = 0
        gap_a = gap_b = False

        for op in wolf:
            if op == MATCH:
                gap_a = gap_b = False
                if i < len(self.__seq1) and j < len(self.__seq2):
                    score += self.__blosum62[self.__seq1[i], self.__seq2[j]] # type:ignore
                    i += 1
                    j += 1

            elif op == GAP_B:
                score += GAP_OPEN if not gap_b else GAP_EXTEND
                gap_b = True
                gap_a = False
                if i < len(self.__seq1):
                    i += 1

            elif op == GAP_A:
                score += GAP_OPEN if not gap_a else GAP_EXTEND
                gap_a = True
                gap_b = False
                if j < len(self.__seq2):
                    j += 1

        while i < len(self.__seq1):
            score += GAP_OPEN if not gap_b else GAP_EXTEND
            gap_b = True
            i += 1

        gap_a = False  # Reset gap state for remaining seq2
        while j < len(self.__seq2):
            score += GAP_OPEN if not gap_a else GAP_EXTEND
            gap_a = True
            j += 1

        return score # type:ignore

    def __update_hierarchy(self) -> None:
        self.__alpha_score = -1e9
        self.__beta_score = -1e9
        self.__delta_score = -1e9

        for i in range(self.__population_size):
            if self.__fitness[i] > self.__alpha_score:
                self.__delta_score = self.__beta_score
                self.__delta = self.__beta.copy()

                self.__beta_score = self.__alpha_score
                self.__beta = self.__alpha.copy()

                self.__alpha_score = self.__fitness[i]
                self.__alpha = self.__wolves[i].copy()
            elif self.__fitness[i] > self.__beta_score:
                self.__delta_score = self.__beta_score
                self.__delta = self.__beta.copy()

                self.__beta_score = self.__fitness[i]
                self.__beta = self.__wolves[i].copy()
            elif self.__fitness[i] > self.__delta_score:
                self.__delta_score = self.__fitness[i]
                self.__delta = self.__wolves[i].copy()

    def __update_positions(self, iteration: int) -> None:
        a = 2.0 - (2.0 * iteration / self.__max_iter)

        for i in range(self.__population_size):
            for j in range(self.__dim):
                r1 = random.random()
                r2 = random.random()

                a1 = 2 * a * r1 - a
                c1 = 2 * r2
                d_alpha = abs(c1 * self.__alpha[j] - self.__wolves[i][j])
                x1 = self.__alpha[j] - a1 * d_alpha

                a2 = 2 * a * random.random() - a
                c2 = 2 * random.random()
                d_beta = abs(c2 * self.__beta[j] - self.__wolves[i][j])
                x2 = self.__beta[j] - a2 * d_beta

                a3 = 2 * a * random.random() - a
                c3 = 2 * random.random()
                d_delta = abs(c3 * self.__delta[j] - self.__wolves[i][j])
                x3 = self.__delta[j] - a3 * d_delta

                new_position: int = round((x1 + x2 + x3) / 3.0)

                self.__wolves[i][j] = max(0, min(2, new_position))
