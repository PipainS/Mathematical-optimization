using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OROM
{
    public enum RestrictionType
    {
        Up,
        Down,
        DoubleRestriction,
        NoneRestriction,
    }

    public class Restriction
    {
        public RestrictionType type { get; set; }
        public double Upper { get; set; }
        public double Lower { get; set; }
        public Restriction(double u, double l, RestrictionType t)
        {
            this.type = t;
            this.Upper = u;
            this.Lower = l;
        }

        public bool isRegion(double x)
        {
            if (type == RestrictionType.Up)
            {
                if (x <= Upper) return true;
                else return false;
            }
            if (type == RestrictionType.Down)
            {
                if (x >= Lower) return true;
                else return false;
            }
            if (type == RestrictionType.DoubleRestriction)
            {
                if (x <= Upper & x >= Lower) return true;
                else return false;
            }
            if (type == RestrictionType.NoneRestriction)
            {
                return true;
            }
            return false;
        }

        public double Proection(double x)
        {
            if (type == RestrictionType.Up)
            {
                if (x > Upper) return Upper;
                else return x;
            }
            if (type == RestrictionType.Down)
            {
                if (x < Lower) return Lower;
                else return x;
            }
            if (type == RestrictionType.DoubleRestriction)
            {
                if (x > Upper) return Upper;
                else if (x < Lower) return Lower;
                else return x;
            }
            if (type == RestrictionType.NoneRestriction)
            {
                return x;
            }
            return x;
        }
    }
    public class RestrictionFunc
    {
        public RestrictionType type { get; set; }
        public double Upper { get; set; }
        public double Lower { get; set; }
        public Func<Vector, double> function;

        public RestrictionFunc(double Upper, double Lower, RestrictionType type, Func<Vector, double> function)
        {
            this.type = type;
            this.Upper = Upper;
            this.Lower = Lower;
            this.function = function;
        }
        public double Proection(double x)
        {
            if (type == RestrictionType.Up)
            {
                if (x > Upper) return Upper;
                else return x;
            }
            if (type == RestrictionType.Down)
            {
                if (x < Lower) return Lower;
                else return x;
            }
            if (type == RestrictionType.DoubleRestriction)
            {
                if (x > Upper) return Upper;
                else if (x < Lower) return Lower;
                else return x;
            }
            if (type == RestrictionType.NoneRestriction)
            {
                return x;
            }
            return x;
        }
        public bool isRegion(Vector x)
        {
            double f = function(x);
            if (type == RestrictionType.Up)
            {
                if (f <= Upper) return true;
                else return false;
            }
            if (type == RestrictionType.Down)
            {
                if (f >= Lower) return true;
                else return false;
            }
            if (type == RestrictionType.DoubleRestriction)
            {
                if (f <= Upper & f >= Lower) return true;
                else return false;
            }
            if (type == RestrictionType.NoneRestriction)
            {
                return true;
            }
            return false;
        }

        public double getValue(Vector x)
        {
            return function(x);
        }

        public double Normalize(Vector v)
        {
            double fx = function(v);
            double norm = fx;
            if (type == RestrictionType.Up)
            {
                if (Upper > 0) norm = fx / Upper;
                if (Upper == 0) norm = fx / Lower + 1.0;
                if (Upper < 0) norm = 2.0 - fx / Upper;
            }
            if (type == RestrictionType.Down)
            {
                if (Lower > 0) norm = 2.0 - fx / Upper;
                if (Upper == 0) norm = 1.0 - fx / Lower;
                if (Upper < 0) norm = fx / Upper;
            }
            if (type == RestrictionType.DoubleRestriction)
            {
                double n1 = (fx - Lower) / (Upper - Lower);
                norm = (Upper - fx) / (Upper - Lower);
                if (n1 > norm) norm = n1;
            }
            return norm;
        }
    }
}
