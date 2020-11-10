package derivative;

import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Date;
import java.util.Map;


public class solver {
    public static void solve(int flag, Map <Character, double[]> args, double precision){
        Date date=new Date();
        long start = date.getTime();
        double [] Aij = args.get('A');
        double [] x = args.get('x');
        double [] n = args.get('n');
        double res = 0;
        if (flag == 1) {
            res = method_chord_f(x, Aij, n, precision);
        } else if (flag == 2) {
            double [] Kij = args.get('K');
            res = method_chord_f2(x, Aij, Kij, n, precision);
        } else if (flag == 3) {
            double [] Kij = args.get('K');
            res = method_chord_f3(x, Aij, Kij, n, precision);
        }
        long finish = date.getTime();
        System.out.println("Exec time: " + (finish - start) + " ms");
        System.out.println("Result: " + res);
    }


    public static double method_chord_f(double [] x, double[] Aij, double [] n, double e) {
        double x_prev = x[0];
        double x_curr = x[1];
        double x_next = 0;
        double n1 = n[0];
        double n2 = n[1];
        double tmp;
        do{
            tmp = x_next;
            x_next = x_curr - f(x_curr, Aij, n1, n2) * (x_prev - x_curr) / (f(x_prev, Aij, n1, n2) - f(x_curr, Aij, n1, n2));
            x_prev = x_curr;
            x_curr = tmp;
        } while (Math.abs(x_next - x_curr) > e);
        return x_next;
    }


    public static double method_chord_f2(double [] x, double[] Aij, double[] Kij, double [] n, double e) {
        double x_prev = x[0];
        double x_curr = x[1];
        double x_next = 0;
        double tmp;
        do{
            tmp = x_next;
            x_next = x_curr - f2(x_curr, Aij, Kij, n) * (x_prev - x_curr) / (f2(x_prev, Aij, Kij, n) - f2(x_curr, Aij, Kij, n));
            x_prev = x_curr;
            x_curr = tmp;
        } while (Math.abs(x_next - x_curr) > e);
        return x_next;
    }

    public static double method_chord_f3(double [] x, double[] Aij, double[] Kij, double [] n, double e) {
        double x_prev = x[0];
        double x_curr = x[1];
        double x_next = 0;
        double tmp;
        do{
            tmp = x_next;
            x_next = x_curr - f3(x_curr, Aij, Kij, n) * (x_prev - x_curr) / (f3(x_prev, Aij, Kij, n) - f3(x_curr, Aij, Kij, n));
            x_prev = x_curr;
            x_curr = tmp;
        } while (Math.abs(x_next - x_curr) > e);
        return x_next;
    }

    private static double f(double x, double @NotNull [] Aij, double n1, double n2){
        double A11 = Aij[0];
        double A12 = Aij[1];
        double A21 = Aij[2];
        double A22 = Aij[3];
        return 1/((n1 - 1)*Math.pow((Math.pow(A11,n1)*A12)/(Math.pow(A12,n1)*(A11 + x)),(1/(n1 - 1)))) -
                1/((n2 - 1  )*Math.pow((Math.pow(A21,n2)*A22)/(Math.pow(A22,n2)*(A21 - x)),(1/(n2 - 1))))+
                1/Math.pow((Math.pow(A11,n1)*A12)/(Math.pow(A12,n1)*(A11 + x)),(1/(n1 - 1)))
                - 1/Math.pow((Math.pow(A21,n2)*A22)/(Math.pow(A22,n2)*(A21 - x)),(1/(n2 - 1)));
    }


    private static double f2(double x, double[] Aij, double[] Kij, double[] ni){
        double A11 = Aij[0];
        double A12 = Aij[1];
        double A21 = Aij[2];
        double A22 = Aij[3];
        double K11 = Kij[0];
        double K12 = Kij[1];
        double K21 = Kij[2];
        double K22 = Kij[3];
        double n1 = ni[0];
        double n2 = ni[1];
        return (A11*A12*K11*K12*n1)/(Math.pow((K11*x + A11*K12),2)*(n1 - 1)*
                Math.pow(((A11*K12)/(K11*x + A11*K12)),((2*n1 - 1)/(n1 - 1)))) -
                (A21*A22*K22*n2*Math.pow((A21 - x),(1/(n2 - 1))))/(K21*(n2 - 1));
    }

    private static double f3(double x, double[] Aij, double[] Kij, double[] ni){
        double A11 = Aij[0];
        double A12 = Aij[1];
        double A21 = Aij[2];
        double A22 = Aij[3];
        double A31 = Aij[4];
        double A32 = Aij[5];
        double K11 = Kij[0];
        double K12 = Kij[1];
        double K21 = Kij[2];
        double K22 = Kij[3];
        double K31 = Kij[4];
        double K32 = Kij[5];
        double n1 = ni[0];
        double n2 = ni[1];
        double n3 = ni[2];
        double v = (A11 * K12) / (K11 * x + A11 * K12);
        double v1 = 1 / Math.pow(v, (n1 / (n1 - 1))) - 1;
        double v2 = A21 * K22 - A12 * K21 * v1;
        return - (A31 * A32 * K32 * n3 * Math.pow((A31 - x),(1/(n3 - 1))))/(K31*(n3 - 1))
                - (A11 * A12 * A21 * A22 * K11 * K12 * K21 * K22 * n1 * n2)/(Math.pow((K11 * x + A11 * K12),2)
                * Math.pow(v2,2) * (n1 - 1) * (n2 - 1) * Math.pow(v,((2 * n1 - 1)/(n1 - 1)))
                * Math.pow(((A21 * K22)/ v2),((2 * n2 - 1)/(n2 - 1))));
    }

}
