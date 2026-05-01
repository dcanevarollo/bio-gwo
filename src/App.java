public class App {

    public static void main(String[] args) {
        final int populationSize = 30;
        final int dimension = 20;
        final int maxIterations = 100;
        final String sequenceA = "ACGTACGT";
        final String sequenceB = "AGTAC";

        GreyWolfOptimizer gwo = new GreyWolfOptimizer(populationSize, dimension, maxIterations, sequenceA, sequenceB);

        gwo.optimize();
        // gwo.printBestAlignment();
    }

}
