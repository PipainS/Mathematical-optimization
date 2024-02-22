using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace OROM
{
    public class Point
    {
        public Vector x { get; set; }
        public double func { get; set; }
        public Point(Vector x, double f)
        {
            this.x = x;
            this.func = f;
        }

        public override string ToString()
        {
            var str = string.Format("{0}: {1}", x, func);
            return str;
        }
    }
}
