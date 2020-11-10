package derivative;

import java.util.HashMap;
import java.util.Map;

public class Main {

    public static void main(String[] args) {
        Map<Character, double[]> params = new HashMap<>();

        params.put('A', new double[]{1232435, 4.3204e9, 1202435, 4.8024e9});
        params.put('x', new double[] {-1, 1});
        params.put('n', new double [] {0.5, 0.6});
        solver.solve(1, params, 0.001);

        params.put('A', new double[]{100, 20000, 1000000, 2000000, 10000, 2000000});
        params.put('x', new double[] {-1, 1});
        params.put('n', new double[] {0.6, 0.5, 0.5});
        params.put('K', new double[] {997, 1000, 997, 1000, 997, 1000});
        solver.solve(2, params, 0.001);
        solver.solve(3, params, 0.001);
    }
}
