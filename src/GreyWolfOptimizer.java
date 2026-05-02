import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class GreyWolfOptimizer {

    private static final int MATCH = 0;

    private static final int GAP_B = 1;

    private static final int GAP_A = 2;

    private static final int GAP_OPEN = -5;

    private static final int GAP_EXTEND = -1;

    private final int populationSize;

    private final int dimension;

    private final int maxIterations;

    private final String sequenceA;

    private final String sequenceB;

    private int[][] wolves;

    private int [] alpha;

    private int [] beta;

    private int [] delta;

    private double[] fitness;

    private double alphaScore = -1e9;

    private double betaScore = -1e9;

    private double deltaScore = -1e9;

    private List<Double> convergence = new ArrayList<>();

    public GreyWolfOptimizer(int populationSize, int dimension, int maxIterations, String sequenceA, String sequenceB) {
        this.populationSize = populationSize;
        this.dimension = dimension;
        this.maxIterations = maxIterations;
        this.sequenceA = sequenceA;
        this.sequenceB = sequenceB;

        this.wolves = new int[populationSize][dimension];
        this.alpha = new int[dimension];
        this.beta = new int[dimension];
        this.delta = new int[dimension];

        this.fitness = new double[populationSize];
    }

    public void optimize() {
        this.initialize();

        for (int iteration = 0; iteration < this.maxIterations; iteration++) {
            for (int i = 0; i < this.populationSize; i++) {
                this.fitness[i] = this.evaluate(this.wolves[i]);
            }

            this.updateHierarchy();
            this.updatePositions(iteration);

            this.convergence.add(this.alphaScore);
        }
    }

    public void printBestAlignment() {
        int[] best = this.alpha;

        StringBuilder alignedA = new StringBuilder();
        StringBuilder alignedB = new StringBuilder();

        int i = 0;
        int j = 0;
        for (int k = 0; k < best.length; k++) {
            int operation = best[k];

            if (operation == MATCH) {
                if (i < this.sequenceA.length() && j < this.sequenceB.length()) {
                    alignedA.append(this.sequenceA.charAt(i++));
                    alignedB.append(this.sequenceB.charAt(j++));
                }
            } else if (operation == GAP_B) {
                if (i < this.sequenceA.length()) {
                    alignedA.append(this.sequenceA.charAt(i++));
                    alignedB.append("-");
                }
            } else if (operation == GAP_A) {
                if (j < this.sequenceB.length()) {
                    alignedA.append("-");
                    alignedB.append(this.sequenceB.charAt(j++));
                }
            }
        }

        System.out.println("\nBest alignment:");
        System.out.println(alignedA);
        System.out.println(alignedB);
        System.out.println("Score: " + this.alphaScore);
    }

    private void initialize() {
        Random random = new Random();

        for (int i = 0; i < populationSize; i++) {
            for (int j = 0; j < dimension; j++) {
                this.wolves[i][j] = random.nextInt(3);
            }
        }
    }

    private int evaluate(int[] wolf) {
        int i = 0;
        int j = 0;
        int score = 0;
        boolean gapA = false;
        boolean gapB = false;

        for (int k = 0; k < wolf.length; k++) {
            int operation = wolf[k];

            if (operation == MATCH) {
                gapA = false;
                gapB = false;

                if (i < this.sequenceA.length() && j < this.sequenceB.length()) {
                    score += Blosum62.score(this.sequenceA.charAt(i++), this.sequenceB.charAt(j++));
                }
            } else if (operation == GAP_B) {
                if (!gapB) {
                    score += GAP_OPEN;
                } else {
                    score += GAP_EXTEND;
                }

                gapB = true;
                gapA = false;

                if (i < this.sequenceA.length()) {
                    i++;
                }
            } else if (operation == GAP_A) {
                if (!gapA) {
                    score += GAP_OPEN;
                } else {
                    score += GAP_EXTEND;
                }

                gapA = true;
                gapB = false;

                if (j < this.sequenceB.length()) {
                    j++;
                }
            }
        }

        return score;
    }

    private void updateHierarchy() {
        this.alphaScore = -1e9;
        this.betaScore = -1e9;
        this.deltaScore = -1e9;

        for (int i = 0; i < this.populationSize; i++) {
            if (this.fitness[i] > this.alphaScore) {
                this.deltaScore = this.betaScore;
                this.delta = this.beta.clone();

                this.betaScore = this.alphaScore;
                this.beta = this.alpha.clone();

                this.alphaScore = this.fitness[i];
                this.alpha = this.wolves[i].clone();
            } else if (this.fitness[i] > this.betaScore) {
                this.deltaScore = this.betaScore;
                this.delta = this.beta.clone();

                this.betaScore = this.fitness[i];
                this.beta = this.wolves[i].clone();
            } else if (this.fitness[i] > this.deltaScore) {
                this.deltaScore = this.fitness[i];
                this.delta = this.wolves[i].clone();
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
                double x1 = this.alpha[j] - a1 * dAlpha;

                double a2 = 2 * a * Math.random() - a;
                double c2 = 2 * Math.random();
                double dBeta = Math.abs(c2 * this.beta[j] - this.wolves[i][j]);
                double x2 = this.beta[j] - a2 * dBeta;

                double a3 = 2 * a * Math.random() - a;
                double c3 = 2 * Math.random();
                double dDelta = Math.abs(c3 * this.delta[j] - this.wolves[i][j]);
                double x3 = this.delta[j] - a3 * dDelta;

                int newPosition = (int) Math.round((x1 + x2 + x3) / 3.0);

                this.wolves[i][j] = Math.max(0, Math.min(2, newPosition));
            }
        }
    }

}
