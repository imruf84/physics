package physics;

public class Numeric {

    public static double[][] cloneArray(double a[][]) {
        double b[][] = new double[a.length][a[0].length];

        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                b[i][j] = a[i][j];
            }
        }

        return b;
    }
}
