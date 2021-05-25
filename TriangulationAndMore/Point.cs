using System;
using System.Collections.Generic;
using System.Text;

namespace TriangulationAndMore
{
    public class Point
    {
        /// <summary>
        /// Used only for generating a unique ID for each instance of this class that gets generated
        /// </summary>
        private static int _counter;

        /// <summary>
        /// Used for identifying an instance of a class; can be useful in troubleshooting when geometry goes weird
        /// (e.g. when trying to identify when Triangle objects are being created with the same Point object twice)
        /// </summary>
        private readonly int _instanceId = _counter++;

        public double X { get; set; }
        public double Y { get; set; }
        public bool IsBoundary { get; set; }
        public double Fi { get; set; }

        public bool IsInnerBoundaryMinus { get; set; }
        public bool IsInnerBoundaryPlus { get; set; }
        public HashSet<Triangle> AdjacentTriangles { get; set; } = new HashSet<Triangle>();

        public Point(double x, double y)
        {
            X = x;
            Y = y;
            IsBoundary = false;
        }

        public Point(double x, double y, double fi)
        {
            X = x;
            Y = y;
            Fi = fi;
            IsBoundary = false;
        }

        public override string ToString()
        {
            // Simple way of seeing what's going on in the debugger when investigating weirdness

            
            return $"{nameof(Point)} {_instanceId} {X:0.##}@{Y:0.##}";
        }

        double sign(Point p1, Point p2, Point p3)
        {
            return (p1.X - p3.X) * (p2.Y - p3.Y) - (p2.X - p3.X) * (p1.Y - p3.Y);
        }

        public bool PointInTriangle(Point pt, Point v1, Point v2, Point v3)
        {
            double d1, d2, d3;
            bool has_neg, has_pos;

            d1 = sign(pt, v1, v2);
            d2 = sign(pt, v2, v3);
            d3 = sign(pt, v3, v1);

            has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            has_pos = (d1 >= 0) || (d2 >= 0) || (d3 >= 0);

            return !(has_neg && has_pos);
        }

    }



    class PointComparer : IEqualityComparer<Point>
    {
        // Products are equal if their names and product numbers are equal.
        public bool Equals(Point x, Point y)
        {

            //Check whether the compared objects reference the same data.
            if (Object.ReferenceEquals(x, y)) return true;

            //Check whether any of the compared objects is null.
            if (Object.ReferenceEquals(x, null) || Object.ReferenceEquals(y, null))
                return false;

            //Check whether the products' properties are equal.
            return x.X == y.X && x.Y == y.Y;
        }

        // If Equals() returns true for a pair of objects
        // then GetHashCode() must return the same value for these objects.

        public int GetHashCode(Point product)
        {
            //Check whether the object is null
            if (Object.ReferenceEquals(product, null)) return 0;

            //Get hash code for the Name field if it is not null.
            int hashX = product.X == null ? 0 : product.X.GetHashCode();
            int hashY = product.Y == null ? 0 : product.Y.GetHashCode();
            //Calculate the hash code for the product.
            return hashX ^ hashY;
        }
    }
}
