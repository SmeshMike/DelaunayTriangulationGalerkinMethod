using System;
using System.Collections.Generic;
using System.Drawing;
using System.Text;
using System.Windows.Media;
using Color = System.Windows.Media.Color;

namespace TriangulationAndMore
{
    public class PointsProcessing
    {
        public class Triangle
        {
            public PointF point1;
            public PointF point2;
            public PointF point3;

            public Triangle(PointF p1, PointF p2, PointF p3)
            {
                point1 = p1;
                point2 = p2;
                point3 = p3;
            }
        }

        private double MaxX { get; set; }
        private double MaxY { get; set; }
        private IEnumerable<Triangle> border;

        public IEnumerable<Point> GeneratePoints(int amount, double maxX, double maxY)
        {
            MaxX = maxX;
            MaxY = maxY;

            // TODO make more beautiful
            var point0 = new Point(0, 0);
            var point1 = new Point(0, MaxY);
            var point2 = new Point(MaxX, MaxY);
            var point3 = new Point(MaxX, 0);
            var points = new List<Point>() { point0, point1, point2, point3 };
            var tri1 = new Triangle(point0, point1, point2);
            var tri2 = new Triangle(point0, point2, point3);
            border = new List<Triangle>() { tri1, tri2 };

            var random = new Random();
            for (int i = 0; i < amount - 4; i++)
            {
                var pointX = random.NextDouble() * MaxX;
                var pointY = random.NextDouble() * MaxY;
                points.Add(new Point(pointX, pointY));
            }

            return points;
        }


        public List<Triangle> triangles;
        public List<PointF> points;

        public List<Triangle> GenerateTriangles(double xMin, double xMax, double yMin, double yMax, int count)
        {
            GenerateDots(xMin, xMax, yMin, yMax, 8, 8);
            FindTriangles();
            return triangles;
        }

        void GenerateDots(double xMin, double xMax, double yMin, double yMax, int width, int height)
        {

            points = new List<PointF>();
            PointF tmpPointF = PointF.Empty;
            float iShift = 0.2f;
            float jShift =0.2f;
            for (int i = -width/2; i < width/2; ++i)
            {
                for (int j = -height/2; j < height/2; ++j)
                {
                    int iSign = i % 2 == 0 ? 1 : -1;
                    int jSign = j % 2 == 0 ? 1 : -1;
                    tmpPointF.X = Convert.ToSingle((i + 1) * 10) + iShift*iSign*jSign + jShift * iSign * jSign;
                    tmpPointF.Y = Convert.ToSingle((j + 1) * 10) + jShift * iSign * jSign +iShift* iSign*jSign;
                    points.Add(tmpPointF);
                }
            }
        }

        void GenerateRandomDots(double xMin, double xMax, double yMin, double yMax, int count)
        {
            var xCoef = Math.Abs(xMin) + Math.Abs(xMax);
            var yCoef = Math.Abs(yMin) + Math.Abs(yMax);

            points = new List<PointF>();
            PointF tmpPointF = PointF.Empty;

            for (int i = 0; i < count; i++)
            {

                tmpPointF.X = (float)(new Random().NextDouble() * xCoef - xCoef / 2);
                tmpPointF.Y = (float)(new Random().NextDouble() * yCoef - yCoef / 2);
                points.Add(tmpPointF);
            }
        }

        void FindTriangles()
        {
            triangles = new List<Triangle>();
            PointF point1, point2, point3 = PointF.Empty;

            for (int i = 0; i < points.Count-2; i++)
            {
                point1 = points[i];
                for (int j = i+1; j < points.Count-1; j++)
                {
                    bool add;
                    point2 = points[j];
                    for (int k = j+1; k < points.Count; k++)
                    {
                        PointF c;
                        float r;
                        FindCircle(points[i], points[j], points[k],out c, out r);
                        point3 = points[k];
                        add = true;

                        for (int l = 0; l < points.Count; l++)
                        {
                            if (l != k && l != j && l != i)
                            {
                                var x = points[l].X;
                                var y = points[l].Y;

                                if (((x - c.X) * (x - c.X) + (y - c.Y) * (y - c.Y) < r * r))
                                {
                                    add = false;
                                    break;
                                }
                            }
                        }
                        if (add)
                            triangles.Add(new Triangle(point1, point2, point3));
                    }
                }
            }
        }



        // Find a circle through the three points.
        private void FindCircle(PointF a, PointF b, PointF c, out PointF center, out float radius)
        {
            // Get the perpendicular bisector of (x1, y1) and (x2, y2).
            var x1 = (b.X + a.X) / 2;
            var y1 = (b.Y + a.Y) / 2;
            var dy1 = b.X - a.X;
            var dx1 = -(b.Y - a.Y);

            // Get the perpendicular bisector of (x2, y2) and (x3, y3).
            var x2 = (c.X + b.X) / 2;
            var y2 = (c.Y + b.Y) / 2;
            var dy2 = c.X - b.X;
            var dx2 = -(c.Y - b.Y);

            // See where the lines intersect.
            FindIntersection(new PointF(x1, y1), new PointF(x1 + dx1, y1 + dy1), new PointF(x2, y2), new PointF(x2 + dx2, y2 + dy2), out var linesIntersect,
                    out var intersection);
            if (!linesIntersect)
            {
                center = new PointF(0, 0);
                radius = 0;
            }
            else
            {
                center = intersection;
                var dx = center.X - a.X;
                var dy = center.Y - a.Y;
                radius = (float) Math.Sqrt(dx * dx + dy * dy);
            }
        }

        private void FindIntersection(PointF p1, PointF p2, PointF p3, PointF p4, out bool linesIntersect, out PointF intersection)
        {
            // Get the segments' parameters.
            var dx12 = p2.X - p1.X;
            var dy12 = p2.Y - p1.Y;
            var dx34 = p4.X - p3.X;
            var dy34 = p4.Y - p3.Y;

            // Solve for t1 and t2
            var denominator = (dy12 * dx34 - dx12 * dy34);

            var t1 = ((p1.X - p3.X) * dy34 + (p3.Y - p1.Y) * dx34) / denominator;
            if (float.IsInfinity(t1))
            {
                // The lines are parallel (or close enough to it).
                linesIntersect = false;
                intersection = new PointF(float.NaN, float.NaN);
                return;
            }

            linesIntersect = true;
            // Find the points of intersection.
            intersection = new PointF(p1.X + dx12 * t1, p1.Y + dy12 * t1);
        }



    }
}
