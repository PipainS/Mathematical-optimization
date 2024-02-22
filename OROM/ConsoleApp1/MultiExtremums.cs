using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OROM
{
    class MultiExtremums
    {
        public static OptimalResult GradientDescent(Vector startPoint, double h, double eps, Func<Vector, double> function)
        {
            var currentPoint = startPoint;
            Vector newPoint = new(currentPoint.Size);
            int iterations = 0;
            while (true)
            {
                iterations++;
                double currentValue = function(currentPoint);
                newPoint = new(currentPoint.Size);
                // производные
                for (var i = 0; i < currentPoint.Size; i++)
                {
                    // Преобразование из многих переменных в одну 
                    Func<double, double> func = x => function(Replace(currentPoint, x, i));
                    newPoint.SetElement(currentPoint[i] - h * (1.0 / startPoint.Size) * GetDerivative(func, currentPoint[i], eps));
                }
                double newValue = function(newPoint);
                // уменьшение шага
                if (newValue > currentValue)
                    h *= 0.8;
                else
                {
                    if (currentValue - newValue <= eps)
                        break;
                    else
                        currentPoint = newPoint;
                }
            }
            double fs = function(newPoint);
            Vector fv = new Vector(new double[] { fs });
            var res = new OptimalResult(newPoint, null, fv, iterations);
            return res;
        }
        private static double GetDerivative(Func<double, double> function, double point, double delta)
        {
            return (function(point + delta) - function(point - delta)) / (2 * delta);
        }

        private static Vector Replace(Vector point, double replace, int replaceIndex)
        {
            var result = new Vector(point.Size);
            for (var i = 0; i < point.Size; i++)
                if (i == replaceIndex)
                    result.SetElement(replace);
                else
                    result.SetElement(point[i]);
            return result;
        }

        public static OptimalResult SteepestDescentMethod(Vector xn, double h, double eps, Func<Vector, double> func)
        {
            int dimension = xn.Size;
            Vector xs = new Vector(dimension);
            int iterationCount = 0;
            double currentObjectiveValue, previousObjectiveValue = 0;
            Vector gradient;
            double delta = eps / 2;

            currentObjectiveValue = func(xn);

            while (Math.Abs(h) > eps)
            {
                iterationCount++;

                // Вычисление градиента
                gradient = new Vector(dimension);
                for (int i = 0; i < dimension; i++)
                {
                    Vector xp = xn.Copy();
                    xp[i] += delta;
                    gradient[i] = (func(xp) - currentObjectiveValue) / delta;
                }

                // Определение оптимального размера шага
                double optimalStepSize = Extremums.StepByStep(0.0, h, eps, h => func(xn - h * gradient));

                // Обновление текущего решения
                xs = xn - (optimalStepSize * gradient);
                previousObjectiveValue = currentObjectiveValue;

                currentObjectiveValue = func(xs);
                Vector displacement = xs - xn;

                // Проверка условия сходимости
                if (displacement.Norma1() < eps)
                {
                    break;
                }

                xn = xs;
            }

            Vector finalValue = new Vector(new double[] { currentObjectiveValue });
            var result = new OptimalResult(xs, null, finalValue, iterationCount);
            return result;
        }


        public static OptimalResult ConjugateGradientMethod(Vector initialVector, double h, double eps, Func<Vector, double> func)
        {
            int dimension = initialVector.Size;
            Vector currentSolution = new Vector(dimension);
            double currentObjectiveValue, previousObjectiveValue = 0, nextObjectiveValue;
            Vector gradient, nextGradient;

            double delta = eps / 2;
            currentObjectiveValue = func(initialVector);

            int iterationCount = 0;
            double weight = 0;
            Vector searchDirection = new Vector(dimension);

            while (Math.Abs(h) > eps)
            {
                iterationCount++;
                
                // Вычисление градиента
                gradient = CalculateGradient(initialVector, delta, func);

                // Выбор направления поиска
                searchDirection = -1 * gradient;

                // Вычисление следующего градиента
                nextGradient = CalculateGradient(initialVector, delta, func);

                // Обновление направления поиска с использованием формулы сопряженных градиентов
                weight = (nextGradient * nextGradient) / (gradient * gradient);
                searchDirection = -1 * nextGradient + weight * searchDirection;

                // Оптимизация шага по направлению поиска
                double optimalStep = Extremums.StepByStep(0.0, h, eps, h => func(initialVector + h * searchDirection));

                // Обновление текущего решения
                currentSolution = initialVector + optimalStep * searchDirection;
                nextObjectiveValue = func(currentSolution);

                // Проверка сходимости
                Vector solutionDifference = currentSolution - initialVector;
                if (solutionDifference.Norma1() < eps)
                    break;

                // Обновление значений для следующей итерации
                initialVector = currentSolution;
                previousObjectiveValue = nextObjectiveValue;
            }

            Vector finalValue = new Vector(new double[] { previousObjectiveValue });
            var result = new OptimalResult(currentSolution, null, finalValue, iterationCount);
            return result;
        }

        private static Vector CalculateGradient(Vector vector, double delta, Func<Vector, double> objectiveFunction)
        {
            int dimension = vector.Size;
            Vector gradient = new Vector(dimension);

            for (int i = 0; i < dimension; i++)
            {
                Vector perturbedVector = vector.Copy();
                perturbedVector[i] += delta;
                gradient[i] = (objectiveFunction(perturbedVector) - objectiveFunction(vector)) / delta;
            }

            return gradient;
        }

        public static OptimalResult RandomMethod(Vector initialVector, double h, double eps, Func<Vector, double> func)
        {
            int dimension = initialVector.Size;
            int numCandidates = 3 * dimension;

            Vector currentSolution = new Vector(dimension);
            double currentObjectiveValue, bestObjectiveValue = 0;
            double[] candidateObjectiveValues = new double[numCandidates];
            Vector[] value = new Vector[numCandidates];

            currentSolution = initialVector;
            currentObjectiveValue = func(initialVector);
            int iterationCount = 0;

            while (Math.Abs(h) > eps)
            {
                iterationCount++;
                double minObjectiveValue = double.MaxValue;
                int bestCandidateIndex = 0;

                for (int candidateIndex = 0; candidateIndex < numCandidates; candidateIndex++)
                {
                    value[candidateIndex] = initialVector + h * Vector.NormalizeRandom(dimension);
                    candidateObjectiveValues[candidateIndex] = func(value[candidateIndex]);

                    if (candidateObjectiveValues[candidateIndex] < minObjectiveValue)
                    {
                        minObjectiveValue = candidateObjectiveValues[candidateIndex];
                        bestCandidateIndex = candidateIndex;
                    }
                }

                currentSolution = value[bestCandidateIndex];
                bestObjectiveValue = func(currentSolution);

                if (bestObjectiveValue > currentObjectiveValue)
                {
                    h = h * 0.5;
                }
                else
                {
                    h = h * 1.2;
                }

                initialVector = currentSolution;
                currentObjectiveValue = bestObjectiveValue;
            }

            Vector finalObjectiveValue = new Vector(new double[] { bestObjectiveValue });
            var result = new OptimalResult(currentSolution, null, finalObjectiveValue, iterationCount);

            return result;
        }


        public static OptimalResult FindMinRandomMethodWithRestrictions(Vector initialPoint, double h, double eps,  Func<Vector, double> func, RestrictionFunc[] restrictionsValue, Restriction[] restrArg)
        {
            int dimension = initialPoint.Size;
            int numTrials = 3 * dimension;

            Vector currentPoint = initialPoint.Copy();
            double currentObjectiveValue, trialObjectiveValue, previousObjectiveValue;
            double[] trialObjectiveValues = new double[numTrials];
            Vector[] trialPoints = new Vector[numTrials];

            currentObjectiveValue = func(initialPoint);
            for (int i = 0; i < dimension; i++)
            {
                if (!restrArg[i].isRegion(currentPoint[i]))
                    currentPoint[i] = restrArg[i].Proection(currentPoint[i]);
            }
            foreach (var restriction in restrictionsValue)
            {
                if (!restriction.isRegion(currentPoint))
                    return null;
            }
            var iterations = 0;
            while (Math.Abs(h) > eps)
            {
                iterations++;
                double minObjectiveValue = double.MaxValue;
                int minIndex = 0;

                for (int i = 0; i < numTrials; i++)
                {
                    do
                    {
                        trialPoints[i] = currentPoint + h * Vector.NormalizeRandom(dimension);

                        bool isValid = restrictionsValue.All(restriction => restriction.isRegion(trialPoints[i]));
                        if (!isValid) continue;

                    } while (false);

                    trialObjectiveValues[i] = func(trialPoints[i]);

                    if (trialObjectiveValues[i] < minObjectiveValue)
                    {
                        minObjectiveValue = trialObjectiveValues[i];
                        minIndex = i;
                    }
                }

                currentPoint = trialPoints[minIndex];
                previousObjectiveValue = currentObjectiveValue;
                currentObjectiveValue = func(currentPoint);

                if (currentObjectiveValue > previousObjectiveValue)
                    h *= 0.5;
                else
                    h *= 1.2;
            }

            Vector constraintValues = new Vector(restrictionsValue.Length);
            for (int i = 0; i < restrictionsValue.Length; i++)
            {
                constraintValues[i] = restrictionsValue[i].getValue(currentPoint);
            }

            double optimalObjectiveValue = func(currentPoint);
            Vector objectiveValues = new Vector(new double[] { optimalObjectiveValue });

            return new OptimalResult(currentPoint, constraintValues, objectiveValues, iterations);
        }


        public static Vector NelderMeadeMethod(Vector xn, double h, double eps, Func<Vector, double> f)
        {
            int n = xn.Size;
            int m = n + 1;
            double alpha = 1;
            double beta = 0.5;
            double gamma = 2.0;
            Vector result = new Vector(n);
            Vector x = xn.Copy();
            Vector[] v = new Vector[m];
            v[0] = x;

            for (int i = 0; i < n; i++)
            {
                Vector xs = x.Copy();
                xs[i] = x[i] + h;
                v[i + 1] = xs;
            }

            double s = getLenghtValue(v);
            while (s > eps)
            {
                Dictionary<Vector, double> functionValues = new Dictionary<Vector, double>();
                foreach (var point in v)
                {
                    double value = f(point);
                    functionValues.Add(point, value);
                }
                var sortedPoints = functionValues.OrderBy(pair => pair.Value).ToDictionary(pair => pair.Key, pair => pair.Value);
                Vector[] sortedKeys = sortedPoints.Keys.ToArray();

                Vector bestPoint = sortedKeys[0];
                Vector secondWorstPoint = sortedKeys[sortedKeys.Length - 2];
                Vector worstPoint = sortedKeys[sortedKeys.Length - 1];

                Vector midpoint = new Vector(n);
                for (int i = 0; i < n; i++)
                    midpoint[i] = (secondWorstPoint[i] + bestPoint[i]) / 2;

                // Отражение
                Vector reflectedPoint = new Vector(n);
                for (int i = 0; i < n; i++)
                    reflectedPoint[i] = midpoint[i] + alpha * (midpoint[i] - worstPoint[i]);

                if (f(reflectedPoint) < f(secondWorstPoint)) { worstPoint = reflectedPoint; }
                else
                {
                    if (f(reflectedPoint) < f(worstPoint)) { worstPoint = reflectedPoint; }
                    Vector contractionPoint = new Vector(n);
                    for (int i = 0; i < n; i++)
                        contractionPoint[i] = (worstPoint[i] + midpoint[i]) / 2;
                    if (f(contractionPoint) < f(worstPoint)) { worstPoint = contractionPoint; }
                }

                // Растяжение
                if (f(reflectedPoint) < f(bestPoint))
                {
                    Vector expansionPoint = new Vector(n);
                    for (int i = 0; i < n; i++)
                        expansionPoint[i] = midpoint[i] + gamma * (reflectedPoint[i] - midpoint[i]);
                    if (f(expansionPoint) < f(reflectedPoint)) { worstPoint = expansionPoint; }
                    else { worstPoint = reflectedPoint; }
                }
                // Сжатие
                if (f(reflectedPoint) < f(secondWorstPoint))
                {
                    Vector contractionPoint = new Vector(n);
                    for (int i = 0; i < n; i++)
                        contractionPoint[i] = midpoint[i] + beta * (worstPoint[i] - midpoint[i]);
                    if (f(contractionPoint) < f(worstPoint)) { worstPoint = contractionPoint; }
                }

                result = bestPoint;
                v[0] = worstPoint; v[v.Length - 2] = secondWorstPoint; v[v.Length - 1] = bestPoint;
                s = getLenghtValue(v);
            }

            return result;
        }

        public static double getLenghtValue(Vector[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double t = (x[i] - x[0]) * (x[i] - x[0]);
                sum += Math.Sqrt(t);
            }
            return sum / x.Length;
        }

        public static OptimalResult Optimize_OZU(Vector initialVector, Restriction[] restrArg, RestrictionFunc[] restrValue, RestrictionFunc[] func, double h, double eps)
        {
            int iterationCount = 0;
            int dimension = initialVector.GetSize();
            int equalityConstraintCount = restrValue.Length;
            int totalConstraints = 3 * dimension;
            Vector currentVector = initialVector.Copy();

            for (int i = 0; i < dimension; i++)
            {
                if (!restrArg[i].isRegion(currentVector[i]))
                    currentVector[i] = restrArg[i].Proection(currentVector[i]);
            }

            bool isFeasible = true;
            for (int i = 0; i < equalityConstraintCount; i++)
            {
                if (!restrValue[i].isRegion(currentVector))
                {
                    isFeasible = false;
                    break;
                }
            }

            if (!isFeasible) return null;

            double currentObjectiveValue = CalculateNormalizedCriterion(currentVector, func);
            Vector bestVector = new Vector(dimension);

            while (Math.Abs(h) > eps)
            {
                double bestObjectiveValue = double.MaxValue;

                for (int i = 0; i < totalConstraints; i++)
                {
                    Vector perturbedVector = new Vector(dimension);
                    bool isPerturbationFeasible = false;

                    while (!isPerturbationFeasible)
                    {
                        perturbedVector = currentVector + h * Vector.NormalizeRandom(dimension);

                        for (int j = 0; j < dimension; j++)
                        {
                            if (!restrArg[j].isRegion(currentVector[j]))
                                perturbedVector[j] = restrArg[j].Proection(perturbedVector[j]);
                        }

                        isPerturbationFeasible = true;

                        for (int j = 0; j < equalityConstraintCount; j++)
                        {
                            if (!restrValue[j].isRegion(perturbedVector))
                            {
                                isPerturbationFeasible = false;
                                break;
                            }
                        }
                    }

                    double perturbedObjectiveValue = CalculateNormalizedCriterion(perturbedVector, func);

                    if (perturbedObjectiveValue < bestObjectiveValue)
                    {
                        bestObjectiveValue = perturbedObjectiveValue;
                        bestVector = perturbedVector.Copy();
                    }
                }

                if (bestObjectiveValue < currentObjectiveValue)
                {
                    currentVector = bestVector.Copy();
                    currentObjectiveValue = bestObjectiveValue;
                    h *= 1.2;
                }
                else
                {
                    h /= 2;
                }

                iterationCount++;
            }

            int objectiveFunctionCount = func.Length;
            Vector optimizedObjectiveValues = new Vector(objectiveFunctionCount);
            for (int i = 0; i < objectiveFunctionCount; i++)
            {
                optimizedObjectiveValues[i] = func[i].getValue(currentVector);
            }

            Vector optimizedEqualityConstraintValues = new Vector(equalityConstraintCount);
            for (int i = 0; i < equalityConstraintCount; i++)
            {
                optimizedEqualityConstraintValues[i] = restrValue[i].getValue(currentVector);
            }

            OptimalResult result = new OptimalResult(currentVector, optimizedEqualityConstraintValues, optimizedObjectiveValues, iterationCount);
            return result;
        }

        public static double CalculateNormalizedCriterion(Vector x, RestrictionFunc[] constraintFunctions)
        {
            double bestCriterion = double.MinValue;

            for (int i = 0; i < constraintFunctions.Length; i++)
            {
                double normalizedCriterion = constraintFunctions[i].Normalize(x);
                if (normalizedCriterion > bestCriterion)
                {
                    bestCriterion = normalizedCriterion;
                }
            }

            return bestCriterion;
        }

    }
}
