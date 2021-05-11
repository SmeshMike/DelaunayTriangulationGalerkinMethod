using System;
using System.Collections.Generic;
using System.Text;
using System.Windows.Media.Media3D;

namespace TriangulationAndMore
{
    public class Triangle3D
    {
        public double A { get; }
        public double B { get; }
        public double C { get; }
        public double D { get; }

        public double S { get; }
        public double S0 { get; }
        public Triangle3D(Point3D m1, Point3D m2, Point3D m3)
        {
            A = GetDetRank2(m2.Y- m1.Y, m3.Y - m1.Y, m2.Z - m1.Z, m3.Z - m1.Z);
            B = -GetDetRank2(m2.X - m1.X, m3.X - m1.X, m2.Z - m1.Z, m3.Z - m1.Z);
            C = GetDetRank2(m2.X - m1.X, m3.X - m1.X, m2.Y - m1.Y, m3.Y - m1.Y);
            D = -m1.X * A - m1.Y * B - m1.Z * C;
            var side1 = Math.Sqrt((m1.X - m2.X) * (m1.X - m2.X) + (m1.Y - m2.Y) * (m1.Y - m2.Y) + (m1.Z - m2.Z) * (m1.Z - m2.Z));
            var side2 = Math.Sqrt((m1.X - m3.X) * (m1.X - m3.X) + (m1.Y - m3.Y) * (m1.Y - m3.Y) + (m1.Z - m3.Z) * (m1.Z - m3.Z));
            var side3 = Math.Sqrt((m1.X - m3.X) * (m1.X - m3.X) + (m1.Y - m3.Y) * (m1.Y - m3.Y) + (m2.Z - m3.Z) * (m2.Z - m3.Z));
            var p = side1 + side2 + side3;
            S = Math.Sqrt(p * (side1 + side2) * (side1 + side3) * (side3 + side2));

            m1.Z = 0;
            m2.Z = 0;
            m3.Z = 0;
            side1 = Math.Sqrt((m1.X - m2.X) * (m1.X - m2.X) + (m1.Y - m2.Y) * (m1.Y - m2.Y) + (m1.Z - m2.Z) * (m1.Z - m2.Z));
            side2 = Math.Sqrt((m1.X - m3.X) * (m1.X - m3.X) + (m1.Y - m3.Y) * (m1.Y - m3.Y) + (m1.Z - m3.Z) * (m1.Z - m3.Z));
            side3 = Math.Sqrt((m1.X - m3.X) * (m1.X - m3.X) + (m1.Y - m3.Y) * (m1.Y - m3.Y) + (m2.Z - m3.Z) * (m2.Z - m3.Z));
            p = side1 + side2 + side3;
            S0 = Math.Sqrt(p * (side1 + side2) * (side1 + side3) * (side3 + side2));
        }

        private double GetDetRank2(double a1, double b1, double a2, double b2)
        {
            return a1 * b2 - b1 * a2;
        }


    }
}
