public class App {

    public static void main(String[] args) {
        final int populationSize = 30;
        final int dimension = 20;
        final int maxIterations = 100;
        final String sequenceA = "HEAGAWGHEE";
        final String sequenceB = "PAWHEAE";

        GreyWolfOptimizer gwo = new GreyWolfOptimizer(populationSize, dimension, maxIterations, sequenceA, sequenceB);

        final long start = System.currentTimeMillis();
        gwo.optimize();
        final long end = System.currentTimeMillis();

        System.out.println("Execution time: " + (end - start) + "ms");

        gwo.printBestAlignment();
    }

}
