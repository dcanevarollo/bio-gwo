public class App {

    public static void main(String[] args) {
        final String sequenceA = "ACGTACGT";
        final String sequenceB = "AGTAC";

        int dimension = sequenceA.length();
        GreyWolfOptimizer gwo = new GreyWolfOptimizer(20, dimension, 100, sequenceA, sequenceB);

        gwo.optimize();
        gwo.printBestAlignment();
    }

}
