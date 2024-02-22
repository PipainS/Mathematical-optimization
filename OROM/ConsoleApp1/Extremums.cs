using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace OROM
{
    public delegate double Function(double x);

    internal class Extremums
    {
        public static double StepByStep(double xn, double h, double eps, Function func)
        {
            double xt = xn;
            double ft = func(xt);
            while (Math.Abs(h) > eps)
            {
                double xs = xt + h;
                double fs = func(xs);
                if (fs < ft)
                {
                    h *= 1.2;
                }
                else
                {
                    h = -h / 2;
                }
                xt = xs;
                ft = fs;
            }
            return xt;
        }

        public static double Gold(double a, double b, double eps, Function f)
        {
            double v = 0.382 * (b - a) + a;
            double w = 0.618 * (b - a) + a;

            double fv = f(v);
            double fw = f(w);
            while (b - a > eps)
            {
                if (fw > fv)
                {
                    b = w;
                    w = v;
                    fw = fv;
                    v = 0.382 * (b - a) + a;
                    fv = f(v);
                }
                else
                {
                    a = v;
                    v = w;
                    fv = fw;
                    w = 0.618 * (b - a) + a;
                    fw = f(w);
                }
            }
            return (a + b) / 2;
        }

        public static double QuadraticApproximation(double xn, double h, double eps, Function f)
        {
            double x1 = xn - h;
            double x2 = xn;
            double x3 = xn + h;

            double y1 = f(x1);
            double y2 = f(x2);
            double y3 = f(x3);

            double xs = double.MinValue;
            double last = double.MaxValue;
            while (last - xs > eps)
            {
                double a = ((y3 - y1) / (x3 - x1) - (y2 - y1) / (x2 - x1)) / (x2 - x1);
                double b = (y2 - y1) / (x2 - x1) - a * (x2 + x1);
                xs = -b / (2 * a);
                double fxs = f(xs);

                if (a > 0)
                {
                    double func = Math.Max(y1, Math.Max(y2, y3));
                    if (func == y1) { last = x1; y1 = fxs; x1 = xs; }
                    else if (func == y2) { last = x2; y2 = fxs; x2 = xs; }
                    else { last = x3; y3 = fxs; x3 = xs; }
                }
                else
                {
                    double func = Math.Min(y1, Math.Min(y2, y3));
                    if (func == y1) { last = x1; y1 = fxs; x1 = xs; }
                    else if (func == y2) { last = y2; y2 = fxs; x2 = xs; }
                    else { last = x3; y3 = fxs; x3 = xs; }
                }
            }
            return xs;
        }

        public static double NewtonsMethod(double currentPoint, double eps, Function func)
        {
            double delta = 0;
            double old_delta = double.MaxValue;
            double deps = eps / 2;
            do
            {
                // Первая производная
                double derivative = (func(currentPoint + deps) - func(currentPoint)) / deps;
                // Вторая производная
                double secondDerivative = (func(currentPoint + deps) - 2 * func(currentPoint) + func(currentPoint - deps)) / (deps * deps);
                // Новое значение
                double newPoint = currentPoint - (derivative / secondDerivative);

                // Абсолютная разность 
                delta = Math.Abs(newPoint - currentPoint);
                currentPoint = newPoint;

                // Проверка на расходимость
                if (delta > old_delta) { return double.NaN; }

                old_delta = delta;
            } while (delta > eps);

            return currentPoint;
        }
    }
}
