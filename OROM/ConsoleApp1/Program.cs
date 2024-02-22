using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;
using System.Numerics;

namespace OROM
{
    class Programm
    {
        public static double FunctionforMnogoug(Vector inp)
        {
            var x = inp[0];
            var y = inp[1];
            return x*x + y*y +2*y +3;
        }
        static void Main(string[] args)
        {
            Console.WriteLine($"Find extremums method step by step:\n" +
                $"equation:  2 * x * x + 5 * x\n" +
                $"Result: " + 
                $"{Extremums.StepByStep(-2, 0.5, 0.001, x => 2 * x * x + 5 * x)}\n");

            Console.WriteLine($"Find extremums method Gold:\n" +
                $"equation:  2 * x * x + 5 * x\n" +
                $"Result: " +
                $"{Extremums.Gold(-2, 0.5, 0.001, x => 2 * x * x + 5 * x)}\n");

            Console.WriteLine($"Find extremums method NewtonsMethod:\n" +
                $"equation:  2 * x * x + 5 * x\n" +
                $"Result: " +
                $"{Extremums.NewtonsMethod(-2, 0.001, x => 2 * x * x + 5 * x)}\n");

            Console.WriteLine($"Find extremums method QuadraticApproximation:\n" +
                $"equation:  Math.Sin(x)\n" +
                $"Result: " +
                $"{Extremums.QuadraticApproximation(-1.5, 0.2, 0.001, x => Math.Sin(x))}\n");

            Vector initialVector = new Vector(new double[] { 3.0, 7.0 }); // Начальная точка поиска
            double stepSize = 1; // Шаг случайного поиска
            double eps = 0.00001; // Критерий останова
            Func<Vector, double> targetFunction = point => point[0] * point[0] +point[1] * point[1] + 2 * point[1] + 3;

            Console.WriteLine("Градиентный метод = {0}\n", MultiExtremums.GradientDescent(initialVector, stepSize, eps, targetFunction));
            Console.WriteLine("Метод случайного поиска = {0}\n", MultiExtremums.RandomMethod(initialVector, stepSize, eps, targetFunction));
            Console.WriteLine("Метод наискорейшего спуска = {0}\n", MultiExtremums.SteepestDescentMethod(initialVector, stepSize, eps, targetFunction));
            Console.WriteLine("Метод сопряженных множеств  = {0}\n", MultiExtremums.ConjugateGradientMethod(initialVector, stepSize, eps, targetFunction));
            Console.WriteLine("Метод многоугольника = {0}\n", MultiExtremums.NelderMeadeMethod(initialVector, stepSize, eps, targetFunction));


            Restriction[] rx = new Restriction[] 
            {
                 new Restriction(3.0, -3.0, RestrictionType.DoubleRestriction),
                 new Restriction(3.0, -3.0, RestrictionType.DoubleRestriction)
            };
            RestrictionFunc[] rf = new RestrictionFunc[] 
            {
                new RestrictionFunc(0.0, 2.0, RestrictionType.Down, x => x[0] + x[1]),

            };
            RestrictionFunc[] rf2 = new RestrictionFunc[]
            {
                new RestrictionFunc(0.0, 5.0, RestrictionType.Up, x => 2 * x[0] + x[1]),
                new RestrictionFunc(2.0, 5.0, RestrictionType.DoubleRestriction, x => x[0] * x[0] + x[1] * x[1])
            };
            Console.WriteLine("Метод ОЗУ = {0}\n", MultiExtremums.Optimize_OZU(initialVector, rx, rf, rf2, stepSize, eps));
        }
    }
}
