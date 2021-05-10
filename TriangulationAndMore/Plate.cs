﻿using System;
using System.Collections.Generic;
using System.Text;
using System.Windows.Media.Media3D;

namespace TriangulationAndMore
{
    public class Plate
    {
        public double A { get; }
        public double B { get; }
        public double C { get; }
        public double D { get; }

        public Plate(Point3D m1, Point3D m2, Point3D m3)
        {
            A = GetDetRank2(m2.Y- m1.Y, m3.Y - m1.Y, m2.Z - m1.Z, m3.Z - m1.Z);
            B = -GetDetRank2(m2.X - m1.X, m3.X - m1.X, m2.Z - m1.Z, m3.Z - m1.Z);
            C = GetDetRank2(m2.X - m1.X, m3.X - m1.X, m2.Y - m1.Y, m3.Y - m1.Y);
            D = -m1.X * A - m1.Y * B - m1.Z * C;
        }

        private double GetDetRank2(double a1, double b1, double a2, double b2)
        {
            return a1 * b2 - b1 * a2;
        }
    }
}
