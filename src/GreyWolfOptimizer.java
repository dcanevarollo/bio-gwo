import java.util.Arrays;

public class GreyWolfOptimizer {

    private final int populationSize;

    private final int dimension;

    private final int maxIterations;

    private final String sequenceA;

    private final String sequenceB;

    private double[][] wolves;

    private double[] fitness;

    private double [] alpha, beta, delta;

    private double alphaScore = Double.NEGATIVE_INFINITY;

    private double betaScore = Double.NEGATIVE_INFINITY;

    private double deltaScore = Double.NEGATIVE_INFINITY;

    public GreyWolfOptimizer(int populationSize, int dimension, int maxIterations, String sequenceA, String sequenceB) {
        this.populationSize = populationSize;
        this.dimension = dimension;
        this.maxIterations = maxIterations;
        this.sequenceA = sequenceA;
        this.sequenceB = sequenceB;

        this.wolves = new double[populationSize][dimension];
        this.fitness = new double[populationSize];

        this.alpha = new double[dimension];
        this.beta = new double[dimension];
        this.delta = new double[dimension];
    }

    public void optimize() {
        this.initialize();

        for (int iteration = 0; iteration < this.maxIterations; iteration++) {
            for (int i = 0; i < this.populationSize; i++) {
                int[] discreteWolf = this.toIntArray(this.wolves[i]);
                this.fitness[i] = this.evaluate(discreteWolf);
            }

            this.updateHierarchy();
            this.updatePositions(iteration);

            if (iteration % 10 == 0) {
                System.out.println("Iteration " + iteration + " | Best score: " + this.alphaScore);
            }
        }
    }

    public void printBestAlignment() {
        int[] best = this.toIntArray(this.alpha);

        StringBuilder alignedA = new StringBuilder();
        StringBuilder alignedB = new StringBuilder();

        int j = 0;
        for (int i = 0; i < best.length; i++) {
            if (best[i] == 1) {
                alignedA.append(this.sequenceA.charAt(i));
                alignedB.append("-");
            } else {
                if (j < this.sequenceB.length()) {
                    alignedA.append(this.sequenceA.charAt(i));
                    alignedB.append(this.sequenceB.charAt(j));
                    j++;
                } else {
                    alignedA.append(this.sequenceA.charAt(i));
                    alignedB.append("-");
                }
            }
        }

        System.out.println("\nBest alignment:");
        System.out.println(alignedA);
        System.out.println(alignedB);
        System.out.println("Score: " + this.alphaScore);
    }

    private void initialize() {
        for (int i = 0; i < populationSize; i++) {
            for (int j = 0; j < dimension; j++) {
                this.wolves[i][j] = Math.random() < 0.5 ? 0 : 1;
            }
        }
    }

    private int[] toIntArray(double[] wolf) {
        int[] array = new int[wolf.length];

        for (int i = 0; i < wolf.length; i++) {
            array[i] = wolf[i] > 0.5 ? 1 : 0;
        }

        return array;
    }

    private int evaluate(int[] wolf) {
        int score = 0;
        int j = 0;

        for (int i = 0; i < wolf.length; i++) {
            if (wolf[i] == 1) {
                score -= 2; // Gap.
            } else {
                if (j >= this.sequenceB.length()) {
                    score -= 2;
                    continue;
                }

                if (this.sequenceA.charAt(i) == this.sequenceB.charAt(j)) {
                    score += 2;
                } else {
                    score -= 1;
                }

                j++;
            }
        }

        return score;
    }

    private void updateHierarchy() {
        for (int i = 0; i < this.populationSize; i++) {
            if (this.fitness[i] > this.alphaScore) {
                this.deltaScore = this.betaScore;
                this.delta = Arrays.copyOf(this.beta, this.dimension);

                this.betaScore = this.alphaScore;
                this.beta = Arrays.copyOf(this.alpha, this.dimension);

                this.alphaScore = this.fitness[i];
                this.alpha = Arrays.copyOf(this.wolves[i], this.dimension);
            } else if (this.fitness[i] > this.betaScore) {
                this.deltaScore = this.betaScore;
                this.delta = Arrays.copyOf(this.beta, this.dimension);

                this.betaScore = this.fitness[i];
                this.beta = Arrays.copyOf(this.wolves[i], this.dimension);
            } else if (this.fitness[i] > this.deltaScore) {
                this.deltaScore = this.fitness[i];
                this.delta = Arrays.copyOf(this.wolves[i], this.dimension);
            }
        }
    }

    private void updatePositions(int iteration) {
        double a = 2.0 - (2.0 - iteration / this.maxIterations);

        for (int i = 0; i < this.populationSize; i++) {
            for (int j = 0; j < this.dimension; j++) {
                double r1 = Math.random();
                double r2 = Math.random();

                double a1 = 2 * a * r1 - a;
                double c1 = 2 * r2;
                double dAlpha = Math.abs(c1 * this.alpha[j] - this.wolves[i][j]);
                double x1 = this.alpha[j] - a1 - dAlpha;

                double a2 = 2 * a * Math.random() - a;
                double c2 = 2 * Math.random();
                double dBeta = Math.abs(c2 * this.beta[j] - this.wolves[i][j]);
                double x2 = this.beta[j] - a2 * dBeta;

                double a3 = 2 * a * Math.random() - a;
                double c3 = 2 * Math.random();
                double dDelta = Math.abs(c3 * this.delta[j] - this.wolves[i][j]);
                double x3 = this.delta[j] - a3 * dDelta;

                double newPosition = (x1 + x2 + x3) / 3.0;

                // Discretization
                this.wolves[i][j] = newPosition > 0.5 ? 1 : 0;
            }
        }
    }

}
